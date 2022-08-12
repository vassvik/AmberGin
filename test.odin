package ambergin

import sdl "vendor:sdl2"

import "core:math"
import "core:fmt"

test :: proc(grid_width, grid_height: i32, velocity_in: [][][2]f32) {
	pressure_ping := make_2D(f32, grid_width, grid_height)
	pressure_pong := make_2D(f32, grid_width, grid_height)
	divergence := make_2D(f32, grid_width, grid_height)
	residual := make_2D(f32, grid_width, grid_height)

	pressure_ping_half := make_2D(f32, grid_width/2, grid_height/2)
	pressure_pong_half := make_2D(f32, grid_width/2, grid_height/2)
	divergence_half := make_2D(f32, grid_width/2, grid_height/2)

	pressure_ping_quarter := make_2D(f32, grid_width/4, grid_height/4)
	pressure_pong_quarter := make_2D(f32, grid_width/4, grid_height/4)
	divergence_quarter := make_2D(f32, grid_width/4, grid_height/4)

	velocity_ping := velocity_in[:]
	velocity_pong := make_2D([2]f32, grid_width, grid_height)

	velocity_encoded := make_2D([2]f32, grid_width, grid_height)
	velocity_encoded_projected := make_2D([2]f32, grid_width, grid_height)
	
	velocity_pong_encoded := make_2D([2]f32, grid_width, grid_height)
	velocity_pong_corrected := make_2D([2]f32, grid_width, grid_height)
	
	velocity_gt_pong := make_2D([2]f32, grid_width, grid_height)
	

	pressure_iterations := 10
	pressure_corrections := 256

	pressure_mode := Pressure_Mode.Vertex
	restriction_mode := Restriction_Mode.Interpolate

	omega := 1.8
	omega_corrections := 1.4
	omega_smooth := 1.001


	iteration_index: int
	iterations: [dynamic]int

	pre_smooth_level0: int = 1
	pre_smooth_level1: int = 1
	solve_iter_level2: int = 2 * int(max(grid_width, grid_height)) / 4
	post_smooth_level1: int = 4
	post_smooth_level0: int = 2


	// Calculate ground truth solution of the raw input
	{
		calc_divergence(velocity_ping, divergence, cast(^Timer(128))nil, pressure_mode)

		// Init Full
		clear_pressure(pressure_ping)

		// Pre-smooth Full
		calc_iterate(&pressure_ping, &pressure_pong, divergence, pre_smooth_level0, 0, omega_smooth, 0.0, pressure_mode)

		for iter in 0..<pressure_iterations {
			// Descend Full -> Quarter
			calc_residual(pressure_ping, pressure_pong, divergence, pressure_mode)
			calc_restriction(pressure_pong, divergence_half, restriction_mode)
			
			// Init Half
			clear_pressure(pressure_ping_half)

			// Pre-smooth Half
			calc_iterate(&pressure_ping_half, &pressure_pong_half, divergence_half, pre_smooth_level1, 0, omega_smooth, 0.0, pressure_mode)

			// Descend Half -> Quarter
			calc_residual(pressure_ping_half, pressure_pong_half, divergence_half, pressure_mode)
			calc_restriction(pressure_pong_half, divergence_quarter, restriction_mode)

			// Init Quarter
			clear_pressure(pressure_ping_quarter)

			// Solve Quarter
			calc_iterate(&pressure_ping_quarter, &pressure_pong_quarter, divergence_quarter, solve_iter_level2, 0, omega, 0.0, pressure_mode)

			// Ascend Quarter -> Half
			calc_prolongation(pressure_ping_quarter, pressure_pong_half)
			for j in 1..<grid_height/2 do for i in 1..<grid_width/2 do pressure_ping_half[j][i] += pressure_pong_half[j][i]

			// Post-smooth Half
			calc_iterate(&pressure_ping_half, &pressure_pong_half, divergence_half, post_smooth_level1, 0, omega_smooth, 0.0, pressure_mode)

			// Ascend Half -> Full
			calc_prolongation(pressure_ping_half, pressure_pong)
			for j in 1..<grid_height do for i in 1..<grid_width do pressure_ping[j][i] += pressure_pong[j][i]

			// Post-smooth Full
			calc_iterate(&pressure_ping, &pressure_pong, divergence, post_smooth_level0, 0, omega_smooth, 0.0, pressure_mode)

			// Corrections
			calc_iterate(&pressure_ping, &pressure_pong, divergence, 0, pressure_corrections, 0.0, omega_corrections, pressure_mode)
		}

		calc_gradient(pressure_ping, velocity_ping, velocity_pong, cast(^Timer(128))nil)

		// Compress and decompress the projected velocity
		for j in 0..<grid_height {
			for i in 0..<grid_width {
				velocity_pong_encoded[j][i] = {
					f32(f16(velocity_pong[j][i].x)),
					f32(f16(velocity_pong[j][i].y)),
				}
			}
		}
	}

	{
		calc_divergence(velocity_pong_encoded, divergence, cast(^Timer(128))nil, pressure_mode)

		// Init Full
		clear_pressure(pressure_ping)

		// Pre-smooth Full
		calc_iterate(&pressure_ping, &pressure_pong, divergence, 4, 0, 4.0/5.0, 0.0, pressure_mode)

		calc_gradient(pressure_ping, velocity_pong_encoded, velocity_pong_corrected, cast(^Timer(128))nil)
	}




	{
		// Calculate solution of the encoded input
		for j in 0..<grid_height {
			for i in 0..<grid_width {
				velocity_encoded[j][i] = {
					f32(f16(velocity_ping[j][i].x)),
					f32(f16(velocity_ping[j][i].y)),
				}
			}
		}

		calc_divergence(velocity_encoded, divergence, cast(^Timer(128))nil, pressure_mode)

		// Init Full
		clear_pressure(pressure_ping)

		// Pre-smooth Full
		calc_iterate(&pressure_ping, &pressure_pong, divergence, pre_smooth_level0, 0, omega_smooth, 0.0, pressure_mode)

		for iter in 0..<pressure_iterations {
			// Descend Full -> Quarter
			calc_residual(pressure_ping, pressure_pong, divergence, pressure_mode)
			calc_restriction(pressure_pong, divergence_half, restriction_mode)
			
			// Init Half
			clear_pressure(pressure_ping_half)

			// Pre-smooth Half
			calc_iterate(&pressure_ping_half, &pressure_pong_half, divergence_half, pre_smooth_level1, 0, omega_smooth, 0.0, pressure_mode)

			// Descend Half -> Quarter
			calc_residual(pressure_ping_half, pressure_pong_half, divergence_half, pressure_mode)
			calc_restriction(pressure_pong_half, divergence_quarter, restriction_mode)

			// Init Quarter
			clear_pressure(pressure_ping_quarter)

			// Solve Quarter
			calc_iterate(&pressure_ping_quarter, &pressure_pong_quarter, divergence_quarter, solve_iter_level2, 0, omega, 0.0, pressure_mode)

			// Ascend Quarter -> Half
			calc_prolongation(pressure_ping_quarter, pressure_pong_half)
			for j in 1..<grid_height/2 do for i in 1..<grid_width/2 do pressure_ping_half[j][i] += pressure_pong_half[j][i]

			// Post-smooth Half
			calc_iterate(&pressure_ping_half, &pressure_pong_half, divergence_half, post_smooth_level1, 0, omega_smooth, 0.0, pressure_mode)

			// Ascend Half -> Full
			calc_prolongation(pressure_ping_half, pressure_pong)
			for j in 1..<grid_height do for i in 1..<grid_width do pressure_ping[j][i] += pressure_pong[j][i]

			// Post-smooth Full
			calc_iterate(&pressure_ping, &pressure_pong, divergence, post_smooth_level0, 0, omega_smooth, 0.0, pressure_mode)

			// Corrections
			calc_iterate(&pressure_ping, &pressure_pong, divergence, 0, pressure_corrections, 0.0, omega_corrections, pressure_mode)
		}

		calc_gradient(pressure_ping, velocity_encoded, velocity_encoded_projected, cast(^Timer(128))nil)
	}

	

	bins_pong_x_abs, bins_pong_y_abs: [32]f64
	bins_pong_x_rel, bins_pong_y_rel: [32]f64
	for j in 0..<grid_height {
		for i in 0..<grid_width {
			v_abs := [2]f32 {
				abs(velocity_pong[j][i].x - velocity_pong_encoded[j][i].x),
				abs(velocity_pong[j][i].y - velocity_pong_encoded[j][i].y),
			}

			v_rel := [2]f32 {
				abs((velocity_pong[j][i].x - velocity_pong_encoded[j][i].x) / velocity_pong[j][i].x),
				abs((velocity_pong[j][i].y - velocity_pong_encoded[j][i].y) / velocity_pong[j][i].y),
			}

			logx_abs := math.log(v_abs.x, 2.0)
			logy_abs := math.log(v_abs.y, 2.0)

			logx_rel := math.log(v_rel.x, 2.0)
			logy_rel := math.log(v_rel.y, 2.0)

			bins_pong_x_abs[clamp(int(logx_abs + 25), 0, 31)] += 1.0
			bins_pong_y_abs[clamp(int(logy_abs + 25), 0, 31)] += 1.0
			bins_pong_x_rel[clamp(int(logx_rel + 25), 0, 31)] += 1.0
			bins_pong_y_rel[clamp(int(logy_rel + 25), 0, 31)] += 1.0
		}
	}

	bins_corr_x_abs, bins_corr_y_abs: [32]f64
	bins_corr_x_rel, bins_corr_y_rel: [32]f64
	for j in 0..<grid_height {
		for i in 0..<grid_width {
			v_abs := [2]f32 {
				abs(velocity_pong[j][i].x - velocity_pong_corrected[j][i].x),
				abs(velocity_pong[j][i].y - velocity_pong_corrected[j][i].y),
			}

			v_rel := [2]f32 {
				abs((velocity_pong[j][i].x - velocity_pong_corrected[j][i].x) / velocity_pong[j][i].x),
				abs((velocity_pong[j][i].y - velocity_pong_corrected[j][i].y) / velocity_pong[j][i].y),
			}

			logx_abs := math.log(v_abs.x, 2.0)
			logy_abs := math.log(v_abs.y, 2.0)

			logx_rel := math.log(v_rel.x, 2.0)
			logy_rel := math.log(v_rel.y, 2.0)

			bins_corr_x_abs[clamp(int(logx_abs + 25), 0, 31)] += 1.0
			bins_corr_y_abs[clamp(int(logy_abs + 25), 0, 31)] += 1.0
			bins_corr_x_rel[clamp(int(logx_rel + 25), 0, 31)] += 1.0
			bins_corr_y_rel[clamp(int(logy_rel + 25), 0, 31)] += 1.0
		}
	}

	bins_encoded_x_abs, bins_encoded_y_abs: [32]f64
	bins_encoded_x_rel, bins_encoded_y_rel: [32]f64
	for j in 0..<grid_height {
		for i in 0..<grid_width {
			v_abs := [2]f32 {
				abs(velocity_pong[j][i].x - velocity_encoded_projected[j][i].x),
				abs(velocity_pong[j][i].y - velocity_encoded_projected[j][i].y),
			}

			v_rel := [2]f32 {
				abs((velocity_pong[j][i].x - velocity_encoded_projected[j][i].x) / velocity_pong[j][i].x),
				abs((velocity_pong[j][i].y - velocity_encoded_projected[j][i].y) / velocity_pong[j][i].y),
			}

			logx_abs := math.log(v_abs.x, 2.0)
			logy_abs := math.log(v_abs.y, 2.0)

			logx_rel := math.log(v_rel.x, 2.0)
			logy_rel := math.log(v_rel.y, 2.0)

			bins_encoded_x_abs[clamp(int(logx_abs + 25), 0, 31)] += 1.0
			bins_encoded_y_abs[clamp(int(logy_abs + 25), 0, 31)] += 1.0
			bins_encoded_x_rel[clamp(int(logx_rel + 25), 0, 31)] += 1.0
			bins_encoded_y_rel[clamp(int(logy_rel + 25), 0, 31)] += 1.0
		}
	}

	for j in 0..<32 {
		fmt.printf("%+ 3d\t% 10f\t% 10f \t% 10f   \t% 10f \t% 10f\t% 10f  \t% 10f \t% 10f \t% 10f  \t% 10f \t% 10f \t% 10f\n", j - 25, 
			       100 * bins_pong_x_abs[j] / f64(grid_width * grid_height), 
			       100 * bins_corr_x_abs[j] / f64(grid_width * grid_height), 
			       100 * bins_encoded_x_abs[j] / f64(grid_width * grid_height), 
			       100 * bins_pong_y_abs[j] / f64(grid_width * grid_height),
			       100 * bins_corr_y_abs[j] / f64(grid_width * grid_height),
			       100 * bins_encoded_y_abs[j] / f64(grid_width * grid_height),
			       100 * bins_pong_x_rel[j] / f64(grid_width * grid_height),
			       100 * bins_corr_x_rel[j] / f64(grid_width * grid_height),
			       100 * bins_encoded_x_rel[j] / f64(grid_width * grid_height),
			       100 * bins_pong_y_rel[j] / f64(grid_width * grid_height),
			       100 * bins_corr_y_rel[j] / f64(grid_width * grid_height),
			       100 * bins_encoded_y_rel[j] / f64(grid_width * grid_height))
	}
}