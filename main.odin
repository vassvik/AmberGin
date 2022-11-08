package ambergin

import sdl "vendor:sdl2"
import sdl_ttf "vendor:sdl2/ttf"
import "core:fmt"
import "core:mem"
import "core:math"
import "core:math/linalg"
import "core:strings"


// Defining these will try to prefer discrete graphics over integrated graphics
@(export, link_name="NvOptimusEnablement")
NvOptimusEnablement: u32 = 0x00000001;

@(export, link_name="AmdPowerXpressRequestHighPerformance")
AmdPowerXpressRequestHighPerformance: i32 = 1;



make_2D :: proc($T: typeid, M, N: i32, allocator := context.allocator) -> [][]T {
	backing := make([]T, int(M)*int(N))
	grid := make([][]T, int(N))
	for row, j in &grid do row = backing[j*int(M):(j+1)*int(M)]
	return grid
}

free_2D :: proc(grid: [][]$T) {
	delete(grid[0])
	delete(grid)
}


main :: proc() {
	grid_width, grid_height: i32 = 255, 255
	// Must be a multiple of 4 minus one, and must be larger than 4
	//  7x7  ->  3x3  -> 1x1
	// 11x11 ->  5x5  -> 2x2
	// 15x15 ->  7x7  -> 3x3
	// 19x19 -> 11x11 -> 4x4
	assert(((grid_width + 1) % 4) == 0)
	assert(((grid_height + 1) % 4) == 0)
	assert(grid_width > 4)
	assert(grid_height > 4)

	// SDL setup
	sdl.Init({.VIDEO, .TIMER})
	defer sdl.Quit()

	sdl_ttf.Init()
	defer sdl_ttf.Quit()

	font := sdl_ttf.OpenFont("C:/Windows/Fonts/Consola.ttf", 16)
	defer sdl_ttf.CloseFont(font)
	assert(font != nil)

	window_width, window_height: i32 = 4*grid_width, 4*grid_height
	window := sdl.CreateWindow("Test SDL", 150, 50, window_width, window_height, {.SHOWN})
	defer sdl.DestroyWindow(window)
	assert(window != nil)

	renderer := sdl.CreateRenderer(window, -1, {.ACCELERATED})
	//renderer := sdl.CreateRenderer(window, -1, {.ACCELERATED, .PRESENTVSYNC})
	defer sdl.DestroyRenderer(renderer)
	assert(renderer != nil)

	display_texture := sdl.CreateTexture(renderer, u32(sdl.PixelFormatEnum.ABGR8888), .STREAMING, 2*grid_width, 2*grid_height)
	defer sdl.DestroyTexture(display_texture)
	assert(display_texture != nil)

	// Grid setup
	density_ping := make_2D(f32, grid_width, grid_height)
	density_pong := make_2D(f32, grid_width, grid_height)

	velocity_x_ping := make_2D(f32, grid_width+1, grid_height+0)
	velocity_x_pong := make_2D(f32, grid_width+1, grid_height+0)
	velocity_y_ping := make_2D(f32, grid_width+0, grid_height+1)
	velocity_y_pong := make_2D(f32, grid_width+0, grid_height+1)

	pressure_ping := make_2D(f32, grid_width+2, grid_height+2)
	pressure_pong := make_2D(f32, grid_width+2, grid_height+2)
	divergence    := make_2D(f32, grid_width+2, grid_height+2)
	residual      := make_2D(f32, grid_width+2, grid_height+2)

	pressure_ping_half := make_2D(f32, (grid_width+1)/2-1+2, (grid_height+1)/2-1+2)
	pressure_pong_half := make_2D(f32, (grid_width+1)/2-1+2, (grid_height+1)/2-1+2)
	divergence_half    := make_2D(f32, (grid_width+1)/2-1+2, (grid_height+1)/2-1+2)

	pressure_ping_quarter := make_2D(f32, (grid_width+1)/4-1+2, (grid_height+1)/4-1+2)
	pressure_pong_quarter := make_2D(f32, (grid_width+1)/4-1+2, (grid_height+1)/4-1+2)
	divergence_quarter    := make_2D(f32, (grid_width+1)/4-1+2, (grid_height+1)/4-1+2)


	// Initialization

	init_sim :: proc(density: [][]f32, velocity_x, velocity_y: [][]f32) {
		grid_height := len(density)
		grid_width := len(density[0])
		
		t := int(10)

		for j in 0..<grid_height {
			for i in 0..<grid_width {
				// Horizontal velocity
				velocity_x[j][i] = 0.0
				velocity_y[j][i] = 0.0
				if j >= grid_height/2-t && j <= grid_height/2+t {
					if i >= grid_width/2-t && i <= grid_height/2+t+1 {
						velocity_x[j][i] = 1.0
					}
				}

				density[j][i] = 0.0
				if j >= grid_height/2-t && j <= grid_height/2+t {
					if i >= grid_width/2-t && i <= grid_width/2+t {
						density[j][i] = 1.0
					}
				}
			}
		}
	}
	init_sim(density_ping, velocity_x_ping, velocity_y_ping)

	omega := 1.86
	omega_smooth := 1.2

	pre_smooth_level0: int = 1
	pre_smooth_level1: int = 1
	solve_iter_level2: int = 1 * int(max(grid_width, grid_height)+1) / 4
	post_smooth_level1: int = 4
	post_smooth_level0: int = 4

	// App state
	should_quit := false
	pressure_iterations := 0

	key_states: map[sdl.Scancode]bool
	key_pressed: map[sdl.Scancode]bool

	frame_timer := create_timer(64)
	start_timer(&frame_timer)

	display_write_timer      := create_timer(64)
	debug_text_timer         := create_timer(64)
	initial_divergence_timer := create_timer(64)
	multigrid_timer          := create_timer(64)
	gradient_timer           := create_timer(64)
	final_divergence_timer   := create_timer(64)


	advect_density :: proc(density_ping, density_pong: ^[][]f32, velocity_x, velocity_y: [][]f32, dt: f32) {
		grid_width := len(density_ping[0])
		grid_height := len(density_ping^)

		for j in 0..<grid_height {
			for i in 0..<grid_width {
				p := [2]f32{f32(i), f32(j)} + 0.5

				vx := 0.5 * (velocity_x[j][i] + velocity_x[j][i+1])
				vy := 0.5 * (velocity_y[j][i] + velocity_y[j+1][i])

				q := p - [2]f32{vx, vy} * dt

				x := math.floor(q.x - 0.5)
				y := math.floor(q.y - 0.5)

				tx := (q.x - 0.5) - x
				ty := (q.y - 0.5) - y

				xi := int(x)
				yi := int(y)

				V00 := 0.0 if xi   < 0 || yi   < 0 || xi   >= grid_width || yi   >= grid_height else density_ping[yi  ][xi  ]
				V10 := 0.0 if xi+1 < 0 || yi   < 0 || xi+1 >= grid_width || yi   >= grid_height else density_ping[yi  ][xi+1]
				V01 := 0.0 if xi   < 0 || yi+1 < 0 || xi   >= grid_width || yi+1 >= grid_height else density_ping[yi+1][xi  ]
				V11 := 0.0 if xi+1 < 0 || yi+1 < 0 || xi+1 >= grid_width || yi+1 >= grid_height else density_ping[yi+1][xi+1]

				V0 := V00 * (1.0 - ty) + V01 * ty
				V1 := V10 * (1.0 - ty) + V11 * ty

				V := V0 * (1.0 - tx) + V1 * tx
				density_pong[j][i] = V
			}
		}

		density_ping^, density_pong^ = density_pong^, density_ping^
	}

	advect_velocity :: proc(velocity_x_ping, velocity_x_pong, velocity_y_ping, velocity_y_pong: ^[][]f32, dt: f32) {
		grid_width := len(velocity_x_ping[0])-1
		grid_height := len(velocity_x_ping^)

		for j in 0..<grid_height {
			for i in 0..<grid_width+1 {

				vx := velocity_x_ping[j][i]
				vy_00 := velocity_y_ping[j  ][i-1] if i > 0 else 0.0
				vy_10 := velocity_y_ping[j  ][i  ] if i < grid_width else 0.0
				vy_01 := velocity_y_ping[j+1][i-1] if i > 0 else 0.0
				vy_11 := velocity_y_ping[j+1][i  ] if i < grid_width else 0.0
				vy := 0.25 * (vy_00 + vy_10 + vy_01 + vy_11)

				p := [2]f32{f32(i), f32(j)} + 0.5
				q := p - [2]f32{vx, vy} * dt

				x := math.floor(q.x - 0.5)
				y := math.floor(q.y - 0.5)

				tx := (q.x - 0.5) - x
				ty := (q.y - 0.5) - y

				xi := int(x)
				yi := int(y)

				V00 := 0.0 if xi   < 0 || yi   < 0 || xi   >= grid_width+1 || yi   >= grid_height else velocity_x_ping[yi  ][xi  ]
				V10 := 0.0 if xi+1 < 0 || yi   < 0 || xi+1 >= grid_width+1 || yi   >= grid_height else velocity_x_ping[yi  ][xi+1]
				V01 := 0.0 if xi   < 0 || yi+1 < 0 || xi   >= grid_width+1 || yi+1 >= grid_height else velocity_x_ping[yi+1][xi  ]
				V11 := 0.0 if xi+1 < 0 || yi+1 < 0 || xi+1 >= grid_width+1 || yi+1 >= grid_height else velocity_x_ping[yi+1][xi+1]

				V0 := V00 * (1.0 - ty) + V01 * ty
				V1 := V10 * (1.0 - ty) + V11 * ty

				V := V0 * (1.0 - tx) + V1 * tx
				velocity_x_pong[j][i] = V
			}
		}

		for j in 0..<grid_height+1 {
			for i in 0..<grid_width {

				vx_00 := velocity_x_ping[j-1][i  ] if j > 0 else 0.0
				vx_10 := velocity_x_ping[j  ][i  ] if j < grid_height else 0.0
				vx_01 := velocity_x_ping[j-1][i+1] if j > 0 else 0.0
				vx_11 := velocity_x_ping[j  ][i+1] if j < grid_height else 0.0
				vx := 0.25 * (vx_00 + vx_10 + vx_01 + vx_11)
				vy := velocity_y_ping[j][i]

				p := [2]f32{f32(i), f32(j)} + 0.5
				q := p - [2]f32{vx, vy} * dt

				x := math.floor(q.x - 0.5)
				y := math.floor(q.y - 0.5)

				tx := (q.x - 0.5) - x
				ty := (q.y - 0.5) - y

				xi := int(x)
				yi := int(y)

				V00 := 0.0 if xi   < 0 || yi   < 0 || xi   >= grid_width || yi   >= grid_height+1 else velocity_y_ping[yi  ][xi  ]
				V10 := 0.0 if xi+1 < 0 || yi   < 0 || xi+1 >= grid_width || yi   >= grid_height+1 else velocity_y_ping[yi  ][xi+1]
				V01 := 0.0 if xi   < 0 || yi+1 < 0 || xi   >= grid_width || yi+1 >= grid_height+1 else velocity_y_ping[yi+1][xi  ]
				V11 := 0.0 if xi+1 < 0 || yi+1 < 0 || xi+1 >= grid_width || yi+1 >= grid_height+1 else velocity_y_ping[yi+1][xi+1]

				V0 := V00 * (1.0 - ty) + V01 * ty
				V1 := V10 * (1.0 - ty) + V11 * ty

				V := V0 * (1.0 - tx) + V1 * tx
				velocity_y_pong[j][i] = V
			}
		}

		velocity_x_ping^, velocity_x_pong^ = velocity_x_pong^, velocity_x_ping^
		velocity_y_ping^, velocity_y_pong^ = velocity_y_pong^, velocity_y_ping^

	}

	for frame := u32(0); !should_quit; frame += 1 {
		stop_timer(&frame_timer)
		start_timer(&frame_timer)

		event: sdl.Event
		clear(&key_pressed)
		for sdl.PollEvent(&event) {
			#partial switch event.type {
			case .QUIT:                             
				should_quit = true
			case .KEYDOWN, .KEYUP:       
				//fmt.println(event.key.keysym)           
				previous_state := key_states[event.key.keysym.scancode]
				if !previous_state && event.key.state == sdl.PRESSED do key_pressed[event.key.keysym.scancode] = true
				key_states[event.key.keysym.scancode] = event.key.state == sdl.PRESSED
			}
		}

		{
			// Input handing

			if key_pressed[.ESCAPE] {
				should_quit = true
			}

			muls: [3]int = {10 if key_states[.LCTRL] else 1, 10 if key_states[.LSHIFT] else 1, 10 if key_states[.LALT] else 1}
			muls2: [3]int = {10 if !key_states[.LCTRL] else 1, 10 if !key_states[.LSHIFT] else 1, 10 if !key_states[.LALT] else 1}
			mul_iter: int = 1
			if key_states[.LCTRL] do mul_iter *= 2
			if key_states[.LSHIFT] do mul_iter *= 2
			if key_states[.LALT] do mul_iter *= -1

			if key_pressed[.NUM1] do pre_smooth_level0  = max(pre_smooth_level0 + mul_iter, 0)
			if key_pressed[.NUM2] do pre_smooth_level1  = max(pre_smooth_level1 + mul_iter, 0)
			if key_pressed[.NUM3] do solve_iter_level2  = max(solve_iter_level2 + mul_iter, 0)
			if key_pressed[.NUM4] do post_smooth_level1 = max(post_smooth_level1 + mul_iter, 0)
			if key_pressed[.NUM5] do post_smooth_level0 = max(post_smooth_level0 + mul_iter, 0)

			if key_pressed[.RIGHT] {
				pressure_iterations += 1*(muls[0])*(muls[1])*(muls[2])
			}
			if key_pressed[.LEFT] {
				pressure_iterations -= 1*(muls[0])*(muls[1])*(muls[2])
			}

			if key_pressed[.K] {
				omega += 1/f64((muls2[0])*(muls2[1])*(muls2[2]))
			}
			if key_pressed[.J] {
				omega -= 1/f64((muls2[0])*(muls2[1])*(muls2[2]))
			}

			if key_pressed[.NUM8] {
				omega_smooth += 1/f64((muls2[0])*(muls2[1])*(muls2[2]))
			}
			if key_pressed[.NUM7] {
				omega_smooth -= 1/f64((muls2[0])*(muls2[1])*(muls2[2]))
			}

			if key_pressed[.R] {
				init_sim(density_ping, velocity_x_ping, velocity_y_ping)
			}

			pressure_iterations = max(0, pressure_iterations)

			omega = clamp(omega, 0.0, 2.0)
		}

		calc_divergence(velocity_x_ping, velocity_y_ping, divergence, &initial_divergence_timer)
		{
			start_timer(&multigrid_timer)
			defer stop_timer(&multigrid_timer)

			// Init Full
			clear_pressure(pressure_ping)

			// Pre-smooth Full
			calc_iterate(&pressure_ping, &pressure_pong, divergence, pre_smooth_level0, omega_smooth)

			for iter in 0..<pressure_iterations {
				// Descend Full -> Quarter
				calc_residual(pressure_ping, pressure_pong, divergence)
				calc_restriction(pressure_pong, divergence_half)
				for j in 0..<(grid_height+1)/2-1 do for i in 0..<(grid_width+1)/2-1 {
					// Free iteration
					pressure_ping_half[1+j][1+i] = f32(omega_smooth) * 0.25 * divergence_half[1+j][1+i]
				}
				
				// Pre-smooth Half
				calc_iterate(&pressure_ping_half, &pressure_pong_half, divergence_half, pre_smooth_level1, omega_smooth)

				// Descend Half -> Quarter
				calc_residual(pressure_ping_half, pressure_pong_half, divergence_half)
				calc_restriction(pressure_pong_half, divergence_quarter)
				for j in 0..<(grid_height+1)/4-1 do for i in 0..<(grid_width+1)/4-1 {
					// Free iteration
					pressure_ping_quarter[1+j][1+i] = f32(omega) * 0.25 * divergence_quarter[1+j][1+i]
				}

				// Init Quarter
				clear_pressure(pressure_ping_quarter)

				// Solve Quarter
				calc_iterate(&pressure_ping_quarter, &pressure_pong_quarter, divergence_quarter, solve_iter_level2, omega)

				// Ascend Quarter -> Half
				calc_prolongation(pressure_ping_quarter, pressure_pong_half)
				for j in 0..<(grid_height+1)/2-1 do for i in 0..<(grid_width+1)/2-1 {
					pressure_ping_half[1+j][1+i] += pressure_pong_half[1+j][1+i]
				}

				// Post-smooth Half
				calc_iterate(&pressure_ping_half, &pressure_pong_half, divergence_half, post_smooth_level1, omega_smooth)

				// Ascend Half -> Full
				calc_prolongation(pressure_ping_half, pressure_pong)
				for j in 0..<grid_height do for i in 0..<grid_width do pressure_ping[1+j][1+i] += pressure_pong[1+j][1+i]

				// Post-smooth Full
				calc_iterate(&pressure_ping, &pressure_pong, divergence, post_smooth_level0, omega_smooth)
			}
		}
		calc_residual(pressure_ping, residual, divergence)
		calc_gradient(pressure_ping, velocity_x_ping, velocity_y_ping, velocity_x_ping, velocity_y_ping, &gradient_timer)
		//calc_divergence(velocity_x_pong, velocity_y_pong, divergence, &final_divergence_timer)

		advect_density(&density_ping, &density_pong, velocity_x_pong, velocity_y_pong, 1.0)
		advect_velocity(&velocity_x_ping, &velocity_x_pong, &velocity_y_ping, &velocity_y_pong, 1.0)
		// Stats
		max_divergence, avg_divergence: f32
		divergence_bins: [9]f64
		for j in 0..<grid_height {
			for i in 0..<grid_width {
				d := abs(divergence[1+j][1+i])

				logd := math.log(d, 10.0)

				max_divergence = max(d, max_divergence)
				avg_divergence += d

				divergence_bins[clamp(int(math.floor(logd)) + 8, 0, 8)] += 1
			}
		}
		avg_divergence /= f32(grid_width)*f32(grid_height)

		if key_pressed[.TAB] {
			sum_bins: f64 = 0
			for i in 0..<9 {
			 	sum_bins += divergence_bins[i]
			}

			fmt.printf("%+ 7f\n", 100*divergence_bins / sum_bins)
		}


		// Stream texture
		{
			start_timer(&display_write_timer)
			defer stop_timer(&display_write_timer)

			// Get texture memory
			locked_pixels_raw: rawptr
			pitch: i32
			sdl.LockTexture(display_texture, nil, &locked_pixels_raw, &pitch)
			defer sdl.UnlockTexture(display_texture)
			assert(locked_pixels_raw != nil)
			assert(pitch == 2*grid_width * size_of([4]u8))

			// Write to pixels
			locked_pixels := mem.slice_ptr(cast(^u32)locked_pixels_raw, int(2*grid_width) * int(2*grid_height))

			velocity_to_u32 :: proc(v: [2]f32) -> u32 {
				r := u32(255 * clamp((0.5 + 0.5*v.x), 0.0, 1.0))
				g := u32(255 * clamp((0.5 + 0.5*v.y), 0.0, 1.0))
				b := u32(0)
				a := u32(255)
				return (r << 0) | (g << 8) | (b << 16) | (a << 24) 
			}

			scalar_to_u32_linear :: proc(scalar: f32) -> u32 {
				r := u32(255 * clamp(scalar, 0.0, 1.0))
				g := u32(255 * clamp(-scalar, 0.0, 1.0))
				b := u32(0)
				a := u32(255)
				return (r << 0) | (g << 8) | (b << 16) | (a << 24) 
			}

			scalar_to_u32_logarithmic :: proc(v: f32) -> u32 {
				logv := math.log(abs(v), 10.0)

				f := math.floor(logv + 7.0)
				i := math.floor(4 * ((logv + 7.0) - f))

				c: [3]f32
			    if (f < 0.0) do c = [3]f32{0.0, 0.0, 0.0};
			    else if (f < 1.0) do c = linalg.mix([3]f32{1.0, 0.0, 0.0}, [3]f32{1.0, 1.0, 1.0}, i / 4.0);
			    else if (f < 2.0) do c = linalg.mix([3]f32{0.0, 1.0, 0.0}, [3]f32{1.0, 1.0, 1.0}, i / 4.0);
			    else if (f < 3.0) do c = linalg.mix([3]f32{0.0, 0.0, 1.0}, [3]f32{1.0, 1.0, 1.0}, i / 4.0);
			    else if (f < 4.0) do c = linalg.mix([3]f32{1.0, 1.0, 0.0}, [3]f32{1.0, 1.0, 1.0}, i / 4.0);
			    else if (f < 5.0) do c = linalg.mix([3]f32{1.0, 0.0, 1.0}, [3]f32{1.0, 1.0, 1.0}, i / 4.0);
			    else if (f < 6.0) do c = linalg.mix([3]f32{0.0, 1.0, 1.0}, [3]f32{1.0, 1.0, 1.0}, i / 4.0);
			    else if (f < 7.0) do c = linalg.mix([3]f32{1.0, 0.5, 0.0}, [3]f32{1.0, 1.0, 1.0}, i / 4.0);
			    else if (f < 8.0) do c = linalg.mix([3]f32{1.0, 1.0, 1.0}, [3]f32{1.0, 1.0, 1.0}, i / 4.0);
			    else do c = [3]f32{1.0, 1.0, 1.0}

				r := u32(255 * clamp(c.x, 0.0, 1.0))
				g := u32(255 * clamp(c.y, 0.0, 1.0))
				b := u32(255 * clamp(c.z, 0.0, 1.0))
				a := u32(255)
				return (r << 0) | (g << 8) | (b << 16) | (a << 24) 
			}


			for y in 0..<u32(grid_height) {
				for x in 0..<u32(grid_width) {
					v := 0.5 * [2]f32{velocity_x_ping[y][x] + velocity_x_ping[y][x+1], velocity_y_ping[y][x] + velocity_y_ping[y+1][x]}
					locked_pixels[u32(2*grid_width)*y + x] = velocity_to_u32(v)
				}
				for x in 0..<u32(grid_width) {
					locked_pixels[u32(2*grid_width)*y + x + u32(grid_width)] = scalar_to_u32_linear(density_ping[y][x])
				}
			}

			for y in 0..<u32(grid_height) {
				for x in 0..<u32(grid_width) {
					locked_pixels[u32(2*grid_width)*(y + u32(grid_height)) + x] = scalar_to_u32_logarithmic(pressure_ping[1+y][1+x])
				}
				for x in 0..<u32(grid_width) {
					locked_pixels[u32(2*grid_width)*(y + u32(grid_height)) + x + u32(grid_width)] = scalar_to_u32_logarithmic(residual[1+y][1+x])
				}
			}
		}

		sdl.RenderClear(renderer);
		sdl.RenderCopy(renderer, display_texture, nil, nil);

		// Render Text	
		if true {
			start_timer(&debug_text_timer)
			defer stop_timer(&debug_text_timer)
			cursor := i32(0)
			text_color := sdl.Color{255, 255, 255, 255};
			render_string(font, renderer, fmt.tprintf("Frame: %d\x00", frame), 0, cursor, text_color); cursor += 24

			render_string(font, renderer, fmt.tprintf("V-Vycles: %d\x00", pressure_iterations), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Iterations: Pre: %d %d, Solve: %d, Post: %d %d\x00", pre_smooth_level0, pre_smooth_level1, solve_iter_level2, post_smooth_level1, post_smooth_level0), 0, cursor, text_color); cursor += 24
			render_string(font, renderer, fmt.tprintf("Omega: Smooth: %f, Solve: %f\x00", omega_smooth, omega), 0, cursor, text_color); cursor += 16

			render_string(font, renderer, fmt.tprintf("Divergence, Max: %.6e, Avg: %.6e\x00", max_divergence, avg_divergence), 0, cursor, text_color); cursor += 24
			
			render_string(font, renderer, fmt.tprintf("Frame              % 7f +/- %f (%f)\x00", frame_timer.average, frame_timer.std, frame_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Display Write      % 7f +/- %f (%f)\x00", display_write_timer.average, display_write_timer.std, display_write_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Initial Divergence % 7f +/- %f (%f)\x00", initial_divergence_timer.average, initial_divergence_timer.std, initial_divergence_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Multigrid          % 7f +/- %f (%f)\x00", multigrid_timer.average, multigrid_timer.std, multigrid_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Gradient           % 7f +/- %f (%f)\x00", gradient_timer.average, gradient_timer.std, gradient_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Final Divergence   % 7f +/- %f (%f)\x00", final_divergence_timer.average, final_divergence_timer.std, final_divergence_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Debug Text         % 7f +/- %f (%f)\x00", debug_text_timer.average, debug_text_timer.std, debug_text_timer.ste), 0, cursor, text_color); cursor += 16
		}
		sdl.RenderPresent(renderer);
	}
}


render_string :: proc(font: ^sdl_ttf.Font, renderer: ^sdl.Renderer, str: string, x, y: i32, color: sdl.Color) {
	// TODO: There seems to be a small memory leak in here? Approximately 100 KB/s
	cstr := strings.unsafe_string_to_cstring(str)

	// TODO: Surely there's an easier way for text that doesn't literally create a surface and require me to make a texture 
	text := sdl_ttf.RenderText_Blended(font, cstr, color)
	defer sdl.FreeSurface(text)
	assert(text != nil)

	text_tex := sdl.CreateTextureFromSurface(renderer, text)
	defer sdl.DestroyTexture(text_tex)
	assert(text_tex != nil)

	rect := sdl.Rect{x, y, text.w, text.h}
	sdl.RenderCopy(renderer, text_tex, nil, &rect);
}

