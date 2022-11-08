package ambergin

import sdl "vendor:sdl2"
import sdl_ttf "vendor:sdl2/ttf"
import "core:fmt"
import "core:mem"
import "core:math"
import "core:math/linalg"
import "core:strings"
import "core:os"


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

	grid_width, grid_height: i32 = 256, 256
	velocity_ping := make_2D([2]f32, grid_width, grid_height)
	velocity_pong := make_2D([2]f32, grid_width, grid_height)
	{
		source_size := min(grid_width, grid_height) / 4
		for j in (grid_height - source_size)/2..<(grid_height + source_size)/2 {
			for i in (grid_width - source_size)/2..<(grid_width + source_size)/2 {
				velocity_ping[j][i].y = 1.0
			}
		}
		/*
		*/
		/*
		for j in 0..<grid_height  {
			for i in 0..<grid_width  {
				x := (f32(i) + 0.5) / f32(grid_width)
				y := (f32(j) + 0.5) / f32(grid_height)
				x *= f32(grid_width) / f32(grid_height)

				w := 0.01 * f32(source_size) / f32(grid_height)
				d := linalg.length([2]f32{x, y} - [2]f32{0.5, 0.5})
				s := 1.0 - math.smoothstep(-w/2, +w/2, d - 1.0/4.0)
				velocity_ping[j][i].x = -2.0*(y - 0.5)*(2.0 * s - 1.0)
				velocity_ping[j][i].y = +2.0*(x - 0.5)*(2.0 * s - 1.0)
			}
		}
		*/
	}

	if false {
		x_velocity_data, x_ok := os.read_entire_file("C:/Users/mv/Downloads/velx_dump")
		y_velocity_data, y_ok := os.read_entire_file("C:/Users/mv/Downloads/vely_dump")
		x_velocity := mem.slice_data_cast([]f32, x_velocity_data)
		y_velocity := mem.slice_data_cast([]f32, y_velocity_data)

		min_x, max_x, min_y, max_y := max(f32), min(f32), max(f32), min(f32)
		for j in 0..<grid_height {
			for i in 0..<grid_width {
				velocity_ping[j][i] = {
					x_velocity[j*grid_width+i] / 256.0,
					y_velocity[j*grid_width+i] / 256.0,
				}
				min_x = min(min_x, x_velocity[j*grid_width+i])
				max_x = max(max_x, x_velocity[j*grid_width+i])
				min_y = min(min_y, y_velocity[j*grid_width+i])
				max_y = max(max_y, y_velocity[j*grid_width+i])
			}
		}
	}

	//test(grid_width, grid_height, velocity_ping)
when true {
	sdl.Init({.VIDEO, .TIMER})
	defer sdl.Quit()

	sdl_ttf.Init()
	defer sdl_ttf.Quit()

	font := sdl_ttf.OpenFont("C:/Windows/Fonts/Consola.ttf", 16)
	defer sdl_ttf.CloseFont(font)
	assert(font != nil)

	window_width, window_height: i32 = 2*512, 2*512
	window := sdl.CreateWindow("Test SDL", 50, 50, window_width, window_height, {.SHOWN})
	defer sdl.DestroyWindow(window)
	assert(window != nil)

	renderer := sdl.CreateRenderer(window, -1, {.ACCELERATED})
	//renderer := sdl.CreateRenderer(window, -1, {.ACCELERATED, .PRESENTVSYNC})
	defer sdl.DestroyRenderer(renderer)
	assert(renderer != nil)


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

	display_texture := sdl.CreateTexture(renderer, u32(sdl.PixelFormatEnum.ABGR8888), .STREAMING, 2*grid_width, 2*grid_height)
	defer sdl.DestroyTexture(display_texture)
	assert(display_texture != nil)


	frame_timer := create_timer(128)
	start_timer(&frame_timer)

	display_write_timer := create_timer(128)
	debug_text_timer := create_timer(128)
	initial_divergence_timer := create_timer(128)
	jacobi_timer := create_timer(128)
	multigrid_timer := create_timer(128)
	gradient_timer := create_timer(128)
	final_divergence_timer := create_timer(128)

	omega := 1.8
	omega_smooth := 1.001

	key_states: map[sdl.Scancode]bool
	key_pressed: map[sdl.Scancode]bool

	Multigrid_Level :: struct {
		level: int,
		iterations: int,
		omega: f64,
	}

	iteration_index: int
	iterations: [dynamic]int

	pre_smooth_level0: int = 1
	pre_smooth_level1: int = 1
	solve_iter_level2: int = 2 * int(max(grid_width, grid_height)) / 4
	post_smooth_level1: int = 4
	post_smooth_level0: int = 2

	should_quit := false
	
	for pressure_iterations := 0; pressure_iterations <= 10; pressure_iterations += 1 do for frame := u32(0); !should_quit; frame += 1 {
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

			pressure_iterations = max(0, pressure_iterations)

			omega = clamp(omega, 0.0, 2.0)
		}

		calc_divergence(velocity_ping, divergence, &initial_divergence_timer)
		if true {
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
				
				// Init Half
				clear_pressure(pressure_ping_half)

				// Pre-smooth Half
				calc_iterate(&pressure_ping_half, &pressure_pong_half, divergence_half, pre_smooth_level1, omega_smooth)

				// Descend Half -> Quarter
				calc_residual(pressure_ping_half, pressure_pong_half, divergence_half)
				calc_restriction(pressure_pong_half, divergence_quarter)

				// Init Quarter
				clear_pressure(pressure_ping_quarter)

				// Solve Quarter
				calc_iterate(&pressure_ping_quarter, &pressure_pong_quarter, divergence_quarter, solve_iter_level2, omega)

				// Ascend Quarter -> Half
				calc_prolongation(pressure_ping_quarter, pressure_pong_half)
				for j in 1..<grid_height/2 do for i in 1..<grid_width/2 do pressure_ping_half[j][i] += pressure_pong_half[j][i]

				// Post-smooth Half
				calc_iterate(&pressure_ping_half, &pressure_pong_half, divergence_half, post_smooth_level1, omega_smooth)

				// Ascend Half -> Full
				calc_prolongation(pressure_ping_half, pressure_pong)
				for j in 1..<grid_height do for i in 1..<grid_width do pressure_ping[j][i] += pressure_pong[j][i]

				// Post-smooth Full
				calc_iterate(&pressure_ping, &pressure_pong, divergence, post_smooth_level0, omega_smooth)
			}


		} else {
			clear_pressure(pressure_ping)
			calc_jacobi(&pressure_ping, &pressure_pong, divergence, &jacobi_timer, pressure_iterations, omega_smooth)
		}
		calc_residual(pressure_ping, residual, divergence)
		calc_gradient(pressure_ping, velocity_ping, velocity_pong, &gradient_timer)
		calc_divergence(velocity_pong, divergence, &final_divergence_timer)

		max_divergence, avg_divergence: f32
		divergence_bins: [9]f64
		for j in 1..<grid_height {
			for i in 1..<grid_width {
				d := abs(divergence[j][i])

				logd := math.log(d, 10.0)

				max_divergence = max(d, max_divergence)
				avg_divergence += d

				divergence_bins[clamp(int(math.floor(logd)) + 8, 0, 8)] += 1
			}
		}
		max_bin, sum_bins: f64 = 0, 0
		for i in 0..<9 {
		 	max_bin = max(max_bin, divergence_bins[i])
		 	sum_bins += divergence_bins[i]
		}

		if key_pressed[.TAB] {
			fmt.printf("%+ 7f\n", 100*divergence_bins / sum_bins)
		}


		precision_bins: [32]f64
		for j in 1..<grid_height {
			for i in 1..<grid_width {
				d := abs(divergence[j][i])
				logd := math.log(d, 2.0)
				precision_bins[clamp(int(logd + 25), 0, 31)] += 1.0
			}
		}

		avg_divergence /= f32(grid_width-1)*f32(grid_height-1)
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
					locked_pixels[u32(2*grid_width)*y + x] = velocity_to_u32(velocity_pong[y][x])
				}
				for x in 0..<u32(grid_width) {
					d := divergence[y][x]
					locked_pixels[u32(2*grid_width)*y + x + u32(grid_width)] = scalar_to_u32_logarithmic(d)
				}
			}

			for y in 0..<u32(grid_height) {
				for x in 0..<u32(grid_width) {
					locked_pixels[u32(2*grid_width)*(y + u32(grid_height)) + x] = scalar_to_u32_logarithmic(pressure_ping[y][x])
				}
				for x in 0..<u32(grid_width) {
					locked_pixels[u32(2*grid_width)*(y + u32(grid_height)) + x + u32(grid_width)] = scalar_to_u32_logarithmic(residual[y][x])
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
			render_string(font, renderer, fmt.tprintf("Frame: %d\x00", frame), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Iterations: %d\x00", pressure_iterations), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Omega: Smooth: %f, Solve: %f\x00", omega_smooth, omega), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Iterations: Pre: %d %d, Solve: %d, Post: %d %d\x00", pre_smooth_level0, pre_smooth_level1, solve_iter_level2, post_smooth_level1, post_smooth_level0), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Frame              % 7f +/- %f (%f)\x00", frame_timer.average, frame_timer.std, frame_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Display Write      % 7f +/- %f (%f)\x00", display_write_timer.average, display_write_timer.std, display_write_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Initial Divergence % 7f +/- %f (%f)\x00", initial_divergence_timer.average, initial_divergence_timer.std, initial_divergence_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Jacobi             % 7f +/- %f (%f)\x00", jacobi_timer.average, jacobi_timer.std, jacobi_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Multigrid          % 7f +/- %f (%f)\x00", multigrid_timer.average, multigrid_timer.std, multigrid_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Gradient           % 7f +/- %f (%f)\x00", gradient_timer.average, gradient_timer.std, gradient_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Final Divergence   % 7f +/- %f (%f)\x00", final_divergence_timer.average, final_divergence_timer.std, final_divergence_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Divergence, Max: %.6e, Avg: %.6e\x00", max_divergence, avg_divergence), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Debug Text         % 7f +/- %f (%f)\x00", debug_text_timer.average, debug_text_timer.std, debug_text_timer.ste), 0, cursor, text_color); cursor += 16
		}
		sdl.RenderPresent(renderer);
	}
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

