package ambergin

import "core:math/linalg"

calc_divergence :: proc(velocity_ping: [][][2]f32, divergence: [][]f32, timer: ^Timer($N)) {
	grid_width := len(velocity_ping[0])
	grid_height := len(velocity_ping)

	start_timer(timer)
	defer stop_timer(timer)

	// Calculate divergence
	for j in 1..<grid_height {
		for i in 1..<grid_width {

			lower_left := velocity_ping[j-1][i-1]
			lower_right := velocity_ping[j-1][i]
			upper_left := velocity_ping[j][i-1]
			upper_right := velocity_ping[j][i]

			vE := 0.5 * (lower_right.x + upper_right.x)
			vW := 0.5 * (lower_left.x + upper_left.x)
			vN := 0.5 * (upper_left.y + upper_right.y)
			vS := 0.5 * (lower_left.y + lower_right.y)

			divergence[j][i] = -((vE - vW) + (vN - vS))
		}
	}
}

clear_pressure :: proc(pressure_ping: [][]f32) {
	grid_width := len(pressure_ping[0])
	grid_height := len(pressure_ping)

	for j in 1..<grid_height {
		for i in 1..<grid_width {
			pressure_ping[j][i] = 0.0
		}
	}
}

calc_stencil_pressure :: proc(input, output, rhs: [][]f32, omega: f64, i, j: int) {
	grid_width := len(input[0])
	grid_height := len(input)

	pE := input[j][i+1] if i < grid_width-1 else 0.0
	pW := input[j][i-1]
	pN := input[j+1][i] if j < grid_height-1 else 0.0
	pS := input[j-1][i]

	p := (rhs[j][i] + pW + pE + pS + pN) / 4.0
	output[j][i] = linalg.mix(input[j][i], p, f32(omega))
}

calc_jacobi :: proc(pressure_ping, pressure_pong: ^[][]f32, divergence: [][]f32, timer: ^Timer($N), iterations: int, omega_solve: f64) {

	start_timer(timer)
	defer stop_timer(timer)

	calc_iterate(pressure_ping, pressure_pong, divergence, iterations, omega_solve)
}

calc_iterate :: proc(pressure_ping, pressure_pong: ^[][]f32, divergence: [][]f32, iterations: int, omega_solve: f64) {
	grid_width := len(pressure_ping[0])
	grid_height := len(pressure_ping)

	// Calculate Jacobi
	for iter := 0; iter < iterations; iter += 1 {
		output := pressure_pong if omega_solve <= 1.0 else pressure_ping
		defer {
			if omega_solve <= 1.0 {
				pressure_ping^, pressure_pong^ = pressure_pong^, pressure_ping^
			}
		}

		{
			for j in 1..<grid_height {
				for i in 1..<grid_width {
					if omega_solve > 1.0 && ((i ~ j) & 1) == 0 do continue
					calc_stencil_pressure(pressure_ping^, output^, divergence, omega_solve, i, j)
				}
			}
			for j in 1..<grid_height {
				for i in 1..<grid_width {
					if omega_solve > 1.0 && ((i ~ j) & 1) == 1 do continue
					calc_stencil_pressure(pressure_ping^, output^, divergence, omega_solve, i, j)
				}
			}
		}
	}
}

calc_residual :: proc(pressure_ping, residual, divergence: [][]f32) {
	grid_width := len(pressure_ping[0])
	grid_height := len(pressure_ping)

	// Calculate Residual
	for j in 1..<grid_height {
		for i in 1..<grid_width {
			pE := pressure_ping[j][i+1] if i < grid_width-1 else 0.0
			pC := pressure_ping[j][i]
			pW := pressure_ping[j][i-1]
			pN := pressure_ping[j+1][i] if j < grid_height-1 else 0.0
			pS := pressure_ping[j-1][i]

			residual[j][i] = divergence[j][i] + pW + pE + pS + pN - 4 * pC
		}
	}
}


calc_gradient :: proc(pressure_ping: [][]f32, velocity_ping, velocity_pong: [][][2]f32, timer: ^Timer($N)) {
	grid_width := len(pressure_ping[0])
	grid_height := len(pressure_ping)

	start_timer(timer)
	defer stop_timer(timer)

	// Calculate gradient
	for j in 0..<grid_height {
		for i in 0..<grid_width {
			lower_left := pressure_ping[j][i]
			lower_right := pressure_ping[j][i+1] if i < grid_width-1 else 0.0
			upper_left := pressure_ping[j+1][i] if j < grid_height-1 else 0.0
			upper_right := pressure_ping[j+1][i+1] if i < grid_width-1 && j < grid_height-1 else 0.0

			pE := 0.5 * (lower_right + upper_right)
			pW := 0.5 * (lower_left + upper_left)
			pN := 0.5 * (upper_left + upper_right)
			pS := 0.5 * (lower_left + lower_right)

			velocity_pong[j][i] = velocity_ping[j][i] - [2]f32{pE - pW, pN - pS}
		}
	}
}


/*
	0 1 2 3 4 5 6 7 8 9 A B C D E F *
	|  \|/ \|/ \|/ \|/ \|/ \|/ \|/  |
	0   1   2   3   4   5   6   7   * Restriction 
	|    \  |  / \  |  / \  |  /    |      |
	|     \ | /   \ | /   \ | /     |      v
	|      \|/     \|/     \|/      | 
	0       1       2       3       *
*/
calc_restriction :: proc(fine_grid: [][]f32, coarse_grid: [][]f32) {
	fine_grid_width := len(fine_grid[0])
	fine_grid_height := len(fine_grid)

	coarse_grid_width := len(coarse_grid[0])
	coarse_grid_height := len(coarse_grid)

	assert(fine_grid_width == coarse_grid_width*2)
	assert(fine_grid_height == coarse_grid_height*2)

	for j in 1..<coarse_grid_height {
		for i in 1..<coarse_grid_width {
			rSW := fine_grid[2*j-1][2*i-1]
			rS  := fine_grid[2*j-1][2*i+0]
			rSE := fine_grid[2*j-1][2*i+1]
			rW  := fine_grid[2*j+0][2*i-1]
			rC  := fine_grid[2*j+0][2*i+0]
			rE  := fine_grid[2*j+0][2*i+1]
			rNW := fine_grid[2*j+1][2*i-1]
			rN  := fine_grid[2*j+1][2*i+0]
			rNE := fine_grid[2*j+1][2*i+1]

			r := (1.0 * (rSW + rSE + rNW + rNE) + 2.0 * (rS + rW + rE + rN) + 4.0 * rC) / 16.0
			coarse_grid[j][i] = 4.0 * r
		}
	}
}

/*
	0       1       2       3       *
	|\     /|\     /|\     /|\     /|
	| \   / | \   / | \   / | \   / | Prolongation
	|  \ /  |  \ /  |  \ /  |  \ /  |      |
	0   1   2   3   4   5   6   7   *      v
	|\ /|\ /|\ /|\ /|\ /|\ /|\ /|\ /|
	0 1 2 3 4 5 6 7 8 9 A B C D E F *
*/
calc_prolongation :: proc(coarse_grid: [][]f32, fine_grid: [][]f32) {
	coarse_grid_width := len(coarse_grid[0])
	coarse_grid_height := len(coarse_grid)
	
	fine_grid_width := len(fine_grid[0])
	fine_grid_height := len(fine_grid)

	assert(fine_grid_width == coarse_grid_width*2)
	assert(fine_grid_height == coarse_grid_height*2)

	for j in 1..<fine_grid_height {
		for i in 1..<fine_grid_width {
			e00 := coarse_grid[(j/2)+0*(j%2)][(i/2)+0*(i%2)]
			e10 := coarse_grid[(j/2)+0*(j%2)][(i/2)+1*(i%2)] if (i/2)+(i%2) < coarse_grid_width else 0.0
			e01 := coarse_grid[(j/2)+1*(j%2)][(i/2)+0*(i%2)] if (j/2)+(j%2) < coarse_grid_height else 0.0
			e11 := coarse_grid[(j/2)+1*(j%2)][(i/2)+1*(i%2)] if (i/2)+(i%2) < coarse_grid_width && (j/2)+(j%2) < coarse_grid_height else 0.0
			fine_grid[j][i] = 0.25 * (e00 + e10 + e01 + e11)
		}
	}
}