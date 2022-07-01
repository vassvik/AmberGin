package ambergin

import "core:math/linalg"

Pressure_Mode :: enum {
	MAC,
	Vertex,
	Wide,
}

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

clear_pressure :: proc(pressure_ping: [][]f32, timer: ^Timer($N)) {
	grid_width := len(pressure_ping[0])
	grid_height := len(pressure_ping)

	start_timer(timer)
	defer stop_timer(timer)

	for j in 1..<grid_height {
		for i in 1..<grid_width {
			pressure_ping[j][i] = 0.0
		}
	}
}

calc_stencil_pressure :: proc(input, output, rhs: [][]f32, omega: f64, pressure_mode: Pressure_Mode, i, j: int) {
	grid_width := len(input[0])
	grid_height := len(input)

	p: f32
	switch pressure_mode {
	case .MAC:
		pE := input[j][i+1] if i < grid_width-1 else 0.0
		pW := input[j][i-1]
		pN := input[j+1][i] if j < grid_height-1 else 0.0
		pS := input[j-1][i]
		p = (rhs[j][i] + pW + pE + pS + pN) / 4.0
	case .Vertex:
		pSW := input[j-1][i-1] 
		pSE := input[j-1][i+1] if i < grid_width-1 else 0.0
		pNW := input[j+1][i-1] if j < grid_height-1 else 0.0
		pNE := input[j+1][i+1] if i < grid_width-1 && j < grid_height-1 else 0.0

		p = (2.0 * rhs[j][i] + pSW + pSE + pNW + pNE) / 4.0
	case .Wide:
		pE  := input[j][i+1] if i < grid_width-1 else 0.0
		pW  := input[j][i-1]
		pN  := input[j+1][i] if j < grid_height-1 else 0.0
		pS  := input[j-1][i]
		pSW := input[j-1][i-1] 
		pSE := input[j-1][i+1] if i < grid_width-1 else 0.0
		pNW := input[j+1][i-1] if j < grid_height-1 else 0.0
		pNE := input[j+1][i+1] if i < grid_width-1 && j < grid_height-1 else 0.0

		p = (3.0 * rhs[j][i] + pSW + pSE + pNW + pNE + pW + pE + pS + pN) / 8.0
	}

	output[j][i] = linalg.mix(input[j][i], p, f32(omega))
}

calc_jacobi :: proc(pressure_ping, pressure_pong: ^[][]f32, divergence: [][]f32, timer: ^Timer($N), 
	                iterations, corrections: int, omega_solve, omega_correction: f64, 
	                pressure_mode: Pressure_Mode) {
	grid_width := len(pressure_ping[0])
	grid_height := len(pressure_ping)

	start_timer(timer)
	defer stop_timer(timer)

	// Calculate Jacobi
	for iter := 0; iter < iterations; iter += 1 {
		output := pressure_pong if omega_solve <= 1.0 else pressure_ping
		defer {
			if omega_solve <= 1.0 {
				pressure_ping^, pressure_pong^ = pressure_pong^, pressure_ping^
			}
		}

		for j in 1..<grid_height {
			for i in 1..<grid_width {
				calc_stencil_pressure(pressure_ping^, output^, divergence, omega_solve, pressure_mode, i, j)
			}
		}
	}

	for iter := 0; iter < corrections; iter += 1 {
		output := pressure_pong if omega_correction <= 1.0 else pressure_ping
		defer if omega_correction <= 1.0 do pressure_ping^, pressure_pong^ = pressure_pong^, pressure_ping^

		for j in 1..<grid_height {
			for i in 1..<grid_width {
				calc_stencil_pressure(pressure_ping^, output^, divergence, omega_correction, Pressure_Mode.Vertex, i, j)
			}
		}
	}
}

calc_residual :: proc(pressure_ping, residual, divergence: [][]f32, timer: ^Timer($N), pressure_mode: Pressure_Mode) {
	grid_width := len(pressure_ping[0])
	grid_height := len(pressure_ping)

	start_timer(timer)
	defer stop_timer(timer)

	// Calculate Residual
	for j in 1..<grid_height {
		for i in 1..<grid_width {
			switch pressure_mode {
			case .MAC:
				pE := pressure_ping[j][i+1] if i < grid_width-1 else 0.0
				pC := pressure_ping[j][i]
				pW := pressure_ping[j][i-1]
				pN := pressure_ping[j+1][i] if j < grid_height-1 else 0.0
				pS := pressure_ping[j-1][i]

				residual[j][i] = divergence[j][i] + pW + pE + pS + pN - 4 * pC
			case .Vertex:
				pSW := pressure_ping[j-1][i-1] 
				pSE := pressure_ping[j-1][i+1] if i < grid_width-1 else 0.0
				pNW := pressure_ping[j+1][i-1] if j < grid_height-1 else 0.0
				pNE := pressure_ping[j+1][i+1] if i < grid_width-1 && j < grid_height-1 else 0.0
				pC  := pressure_ping[j][i]

				residual[j][i] = 2.0 * divergence[j][i] + pSW + pSE + pNW + pNE - 4 * pC
			case .Wide:
				pE  := pressure_ping[j][i+1] if i < grid_width-1 else 0.0
				pW  := pressure_ping[j][i-1]
				pN  := pressure_ping[j+1][i] if j < grid_height-1 else 0.0
				pS  := pressure_ping[j-1][i]
				pSW := pressure_ping[j-1][i-1] 
				pSE := pressure_ping[j-1][i+1] if i < grid_width-1 else 0.0
				pNW := pressure_ping[j+1][i-1] if j < grid_height-1 else 0.0
				pNE := pressure_ping[j+1][i+1] if i < grid_width-1 && j < grid_height-1 else 0.0
				pC  := pressure_ping[j][i]

				residual[j][i] = 3.0 * divergence[j][i] + pSW + pSE + pNW + pNE + pW + pE + pS + pN - 8.0 * pC
			}
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
