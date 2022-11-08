package ambergin

import "core:math/linalg"

calc_divergence :: proc(velocity_x, velocity_y: [][]f32, divergence: [][]f32) {
	grid_width := len(velocity_x[0]) - 1
	grid_height := len(velocity_x)
	
	// Calculate divergence
	for j in 0..<grid_height {
		for i in 0..<grid_width {
			left  := velocity_x[j  ][i  ]
			right := velocity_x[j  ][i+1]
			down  := velocity_y[j  ][i  ]
			up    := velocity_y[j+1][i  ]

			divergence[1+j][1+i] = -((right - left) + (up - down))
		}
	}
}

clear_pressure :: proc(pressure_ping: [][]f32) {
	grid_width := len(pressure_ping[0]) - 2
	grid_height := len(pressure_ping) - 2

	for j in 0..<grid_height {
		for i in 0..<grid_width {
			pressure_ping[1+j][1+i] = 0.0
		}
	}
}

calc_stencil_pressure :: proc(input, output, rhs: [][]f32, omega: f64, i, j: int) {
	pE := input[j][i+1]
	pW := input[j][i-1]
	pN := input[j+1][i]
	pS := input[j-1][i]

	p := (rhs[j][i] + pW + pE + pS + pN) / 4.0
	output[j][i] = linalg.mix(input[j][i], p, f32(omega))
}

calc_iterate :: proc(pressure_ping, pressure_pong: ^[][]f32, divergence: [][]f32, iterations: int, omega_solve: f64) {
	grid_width := len(pressure_ping[0]) - 2
	grid_height := len(pressure_ping) - 2

	// Calculate Jacobi
	for iter := 0; iter < iterations; iter += 1 {
		output := pressure_pong if omega_solve <= 1.0 else pressure_ping
		defer {
			if omega_solve <= 1.0 {
				pressure_ping^, pressure_pong^ = pressure_pong^, pressure_ping^
			}
		}

		if omega_solve > 1.0 {
			for j in 0..<grid_height {
				for i in 0..<grid_width {
					if ((i ~ j) & 1) == 0 do continue
					calc_stencil_pressure(pressure_ping^, output^, divergence, omega_solve, 1+i, 1+j)
				}
			}
			for j in 0..<grid_height {
				for i in 0..<grid_width {
					if ((i ~ j) & 1) == 1 do continue
					calc_stencil_pressure(pressure_ping^, output^, divergence, omega_solve, 1+i, 1+j)
				}
			}
		} else {
			for j in 0..<grid_height {
				for i in 0..<grid_width {
					calc_stencil_pressure(pressure_ping^, output^, divergence, omega_solve, 1+i, 1+j)
				}
			}
		}
	}
}

calc_residual :: proc(pressure_ping, residual, divergence: [][]f32) {
	grid_width := len(pressure_ping[0]) - 2
	grid_height := len(pressure_ping) - 2

	// Calculate Residual
	for j in 0..<grid_height {
		for i in 0..<grid_width {
			pE := pressure_ping[1+j  ][1+i+1]
			pC := pressure_ping[1+j  ][1+i  ]
			pW := pressure_ping[1+j  ][1+i-1]
			pN := pressure_ping[1+j+1][1+i  ]
			pS := pressure_ping[1+j-1][1+i  ]

			residual[1+j][1+i] = divergence[1+j][1+i] + pW + pE + pS + pN - 4 * pC
		}
	}
}


calc_gradient :: proc(pressure_ping, velocity_x_ping, velocity_y_ping, velocity_x_pong, velocity_y_pong: [][]f32) {
	grid_width := len(pressure_ping[0]) - 2
	grid_height := len(pressure_ping) - 2

	// Calculate gradient x
	for j in 0..<grid_height {
		for i in 0..<grid_width+1 {
			left  := pressure_ping[1+j][i ]
			right := pressure_ping[1+j][i+1]

			velocity_x_pong[j][i] = velocity_x_ping[j][i] - (right - left)
		}
	}

	for j in 0..<grid_height+1 {
		for i in 0..<grid_width {
			down := pressure_ping[j  ][1+i]
			up   := pressure_ping[j+1][1+i]

			velocity_y_pong[j][i] = velocity_y_ping[j][i] - (up - down)
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
	coarse_grid_width := len(coarse_grid[0]) - 2
	coarse_grid_height := len(coarse_grid) - 2

	for j in 1..=coarse_grid_height {
		for i in 1..=coarse_grid_width {
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
	fine_grid_width := len(fine_grid[0]) - 2
	fine_grid_height := len(fine_grid) - 2

	for j in 1..=fine_grid_height {
		for i in 1..=fine_grid_width {
			e00 := coarse_grid[(j/2)+0*(j%2)][(i/2)+0*(i%2)]
			e10 := coarse_grid[(j/2)+0*(j%2)][(i/2)+1*(i%2)]
			e01 := coarse_grid[(j/2)+1*(j%2)][(i/2)+0*(i%2)]
			e11 := coarse_grid[(j/2)+1*(j%2)][(i/2)+1*(i%2)]
			fine_grid[j][i] = 0.25 * (e00 + e10 + e01 + e11)
		}
	}
}