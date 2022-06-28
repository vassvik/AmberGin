package ambergin

import sdl "vendor:sdl2"
import sdl_ttf "vendor:sdl2/ttf"
import "core:fmt"
import "core:mem"
import "core:math"
import "core:math/linalg"
import "core:strings"

main :: proc() {
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

	grid_width, grid_height: i32 = 32, 32
	velocity_ping := make_2D([2]f32, grid_width, grid_height)
	velocity_pong := make_2D([2]f32, grid_width, grid_height)
	{
		source_size := min(grid_width, grid_height) / 4
		for j in (grid_height - source_size)/2..<(grid_height + source_size)/2 {
			for i in (grid_width - source_size)/2..<(grid_width + source_size)/2 {
				velocity_ping[j][i].y = 1.0
			}
		}
	}

	divergence := make_2D(f32, grid_width, grid_height)
	residual := make_2D(f32, grid_width, grid_height)
	pressure_ping := make_2D(f32, grid_width, grid_height)
	pressure_pong := make_2D(f32, grid_width, grid_height)

	display_texture := sdl.CreateTexture(renderer, u32(sdl.PixelFormatEnum.ABGR8888), .STREAMING, 2*grid_width, 2*grid_height)
	defer sdl.DestroyTexture(display_texture)
	assert(display_texture != nil)


	frame_timer := create_timer(128)
	start_timer(&frame_timer)

	display_write_timer := create_timer(128)
	debug_text_timer := create_timer(128)
	initial_divergence_timer := create_timer(128)
	clear_pressure_timer := create_timer(128)
	jacobi_timer := create_timer(128)
	gradient_timer := create_timer(128)
	residual_timer := create_timer(128)
	final_divergence_timer := create_timer(128)

	pressure_iterations := 1
	pressure_corrections := 0

	Pressure_Mode :: enum {
		MAC,
		Vertex,
		Wide,
	}
	pressure_mode := Pressure_Mode.MAC

	omega := 1.8
	omega_corrections := 0.8

	key_states: map[sdl.Scancode]bool
	key_pressed: map[sdl.Scancode]bool



	should_quit := false
	for frame := u32(0); !should_quit; frame += 1 {
		stop_timer(&frame_timer)
		start_timer(&frame_timer)

		event: sdl.Event
		clear(&key_pressed)
		for sdl.PollEvent(&event) != 0 {
			#partial switch event.type {
			case .QUIT:                             
				should_quit = true
			case .KEYDOWN, .KEYUP:       
			fmt.println(event.key.keysym)           
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

			if key_pressed[.RIGHT] {
				pressure_iterations += 1*(muls[0])*(muls[1])*(muls[2])
			}
			if key_pressed[.LEFT] {
				pressure_iterations -= 1*(muls[0])*(muls[1])*(muls[2])
			}

			if key_pressed[.UP] {
				pressure_corrections += 1*(muls[0])*(muls[1])*(muls[2])
			}
			if key_pressed[.DOWN] {
				pressure_corrections -= 1*(muls[0])*(muls[1])*(muls[2])
			}

			if key_pressed[.K] {
				omega += 1/f64((muls2[0])*(muls2[1])*(muls2[2]))
			}
			if key_pressed[.J] {
				omega -= 1/f64((muls2[0])*(muls2[1])*(muls2[2]))
			}
			if key_pressed[.O] {
				omega_corrections += 1/f64((muls2[0])*(muls2[1])*(muls2[2]))
			}
			if key_pressed[.I] {
				omega_corrections -= 1/f64((muls2[0])*(muls2[1])*(muls2[2]))
			}

			if key_pressed[.M] {
				pressure_mode = Pressure_Mode((int(pressure_mode) + 1) % len(Pressure_Mode))
			}

			pressure_iterations = max(0, pressure_iterations)
			pressure_corrections = max(0, pressure_corrections)

			omega = clamp(omega, 0.0, 2.0)
			omega_corrections = clamp(omega_corrections, 0.0, 2.0)
		}

		{
			start_timer(&initial_divergence_timer)
			defer stop_timer(&initial_divergence_timer)

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
		{
			start_timer(&clear_pressure_timer)
			defer stop_timer(&clear_pressure_timer)

			for j in 1..<grid_height {
				for i in 1..<grid_width {
					pressure_ping[j][i] = 0.0
				}
			}
		}
		{
			start_timer(&jacobi_timer)
			defer stop_timer(&jacobi_timer)

			// Calculate Jacobi
			for iter := 0; iter < pressure_iterations; iter += 1{
				output := pressure_pong if omega <= 1.0 else pressure_ping
				defer {
					if omega <= 1.0 {
						pressure_ping, pressure_pong = pressure_pong, pressure_ping
					}
				}

				for j in 1..<grid_height {
					for i in 1..<grid_width {
						switch pressure_mode {
						case .MAC:
							pE := pressure_ping[j][i+1] if i < grid_width-1 else 0.0
							pW := pressure_ping[j][i-1]
							pN := pressure_ping[j+1][i] if j < grid_height-1 else 0.0
							pS := pressure_ping[j-1][i]
							pC  := pressure_ping[j][i]

							p := (divergence[j][i] + pW + pE + pS + pN) / 4.0
							output[j][i] = linalg.mix(pC, p, f32(omega))
						case .Vertex:
							pSW := pressure_ping[j-1][i-1] 
							pSE := pressure_ping[j-1][i+1] if i < grid_width-1 else 0.0
							pNW := pressure_ping[j+1][i-1] if j < grid_height-1 else 0.0
							pNE := pressure_ping[j+1][i+1] if i < grid_width-1 && j < grid_height-1 else 0.0
							pC  := pressure_ping[j][i]

							p := (2.0 * divergence[j][i] + pSW + pSE + pNW + pNE) / 4.0
							output[j][i] = linalg.mix(pC, p, f32(omega))
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

							p := (3.0 * divergence[j][i] + pSW + pSE + pNW + pNE + pW + pE + pS + pN) / 8.0
							output[j][i] = linalg.mix(pC, p, f32(omega))
						}
					}
				}
			}

			for iter := 0; iter < pressure_corrections; iter += 1 {
				output := pressure_pong if omega_corrections <= 1.0 else pressure_ping
				defer if omega_corrections <= 1.0 do pressure_ping, pressure_pong = pressure_pong, pressure_ping

				for j in 1..<grid_height {
					for i in 1..<grid_width {
						switch Pressure_Mode.Vertex {
						case .MAC:
							pE := pressure_ping[j][i+1] if i < grid_width-1 else 0.0
							pW := pressure_ping[j][i-1]
							pN := pressure_ping[j+1][i] if j < grid_height-1 else 0.0
							pS := pressure_ping[j-1][i]
							pC  := pressure_ping[j][i]

							p := (divergence[j][i] + pW + pE + pS + pN) / 4.0
							output[j][i] = linalg.mix(pC, p, f32(omega_corrections))
						case .Vertex:
							pSW := pressure_ping[j-1][i-1] 
							pSE := pressure_ping[j-1][i+1] if i < grid_width-1 else 0.0
							pNW := pressure_ping[j+1][i-1] if j < grid_height-1 else 0.0
							pNE := pressure_ping[j+1][i+1] if i < grid_width-1 && j < grid_height-1 else 0.0
							pC  := pressure_ping[j][i]

							p := (2.0 * divergence[j][i] + pSW + pSE + pNW + pNE) / 4.0
							output[j][i] = linalg.mix(pC, p, f32(omega_corrections))
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

							p := (3.0 * divergence[j][i] + pSW + pSE + pNW + pNE + pW + pE + pS + pN) / 8.0
							output[j][i] = linalg.mix(pC, p, f32(omega_corrections))
						}
					}
				}
			}
		}

		{
			start_timer(&residual_timer)
			defer stop_timer(&residual_timer)

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

						r := divergence[j][i] + pW + pE + pS + pN - 4 * pC
						residual[j][i] = r
					case .Vertex:
						pSW := pressure_ping[j-1][i-1] 
						pSE := pressure_ping[j-1][i+1] if i < grid_width-1 else 0.0
						pNW := pressure_ping[j+1][i-1] if j < grid_height-1 else 0.0
						pNE := pressure_ping[j+1][i+1] if i < grid_width-1 && j < grid_height-1 else 0.0
						pC  := pressure_ping[j][i]

						r := 2.0 * divergence[j][i] + pSW + pSE + pNW + pNE - 4 * pC
						residual[j][i] = r
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

						r := 3.0 * divergence[j][i] + pSW + pSE + pNW + pNE + pW + pE + pS + pN - 8.0 * pC
						residual[j][i] = r
					}
				}
			}
		}

		{
			start_timer(&gradient_timer)
			defer stop_timer(&gradient_timer)

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

		{
			start_timer(&final_divergence_timer)
			defer stop_timer(&final_divergence_timer)

			// Calculate divergence
			for j in 1..<grid_height {
				for i in 1..<grid_width {

					lower_left := velocity_pong[j-1][i-1]
					lower_right := velocity_pong[j-1][i]
					upper_left := velocity_pong[j][i-1]
					upper_right := velocity_pong[j][i]

					vE := 0.5 * (lower_right.x + upper_right.x)
					vW := 0.5 * (lower_left.x + upper_left.x)
					vN := 0.5 * (upper_left.y + upper_right.y)
					vS := 0.5 * (lower_left.y + lower_right.y)

					divergence[j][i] = -((vE - vW) + (vN - vS))
				}
			}
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
				r := u32(255 * (0.5 + 0.5*v.x))
				g := u32(255 * (0.5 + 0.5*v.y))
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
					locked_pixels[u32(2*grid_width)*y + x + u32(grid_width)] = scalar_to_u32_logarithmic(divergence[y][x])
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
		{
			start_timer(&debug_text_timer)
			defer stop_timer(&debug_text_timer)
			cursor := i32(0)
			text_color := sdl.Color{255, 255, 255, 255};
			render_string(font, renderer, fmt.tprintf("Frame: %d\x00", frame), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Iterations: %d (%d), mode: %v, omega: %f (%f)\x00", pressure_iterations, pressure_corrections, pressure_mode, omega, omega_corrections), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Frame              % 7f +/- %f (%f)\x00", frame_timer.average, frame_timer.std, frame_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Display Write      % 7f +/- %f (%f)\x00", display_write_timer.average, display_write_timer.std, display_write_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Initial Divergence % 7f +/- %f (%f)\x00", initial_divergence_timer.average, initial_divergence_timer.std, initial_divergence_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Clear Puressure    % 7f +/- %f (%f)\x00", clear_pressure_timer.average, clear_pressure_timer.std, clear_pressure_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Jacobi             % 7f +/- %f (%f)\x00", jacobi_timer.average, jacobi_timer.std, jacobi_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Residual           % 7f +/- %f (%f)\x00", residual_timer.average, residual_timer.std, residual_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Gradient           % 7f +/- %f (%f)\x00", gradient_timer.average, gradient_timer.std, gradient_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Final Divergence   % 7f +/- %f (%f)\x00", final_divergence_timer.average, final_divergence_timer.std, final_divergence_timer.ste), 0, cursor, text_color); cursor += 16
			render_string(font, renderer, fmt.tprintf("Debug Text         % 7f +/- %f (%f)\x00", debug_text_timer.average, debug_text_timer.std, debug_text_timer.ste), 0, cursor, text_color); cursor += 16
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

