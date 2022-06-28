package ambergin

import sdl "vendor:sdl2"
import sdl_ttf "vendor:sdl2/ttf"
import "core:fmt"
import "core:mem"
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

	should_quit := false
	for frame := u32(0); !should_quit; frame += 1 {
		stop_timer(&frame_timer)
		start_timer(&frame_timer)

		event: sdl.Event
		for sdl.PollEvent(&event) != 0 {
			#partial switch event.type {
			case .QUIT:                             
				should_quit = true
			case .KEYDOWN, .KEYUP:                  
				//fmt.println(i, event.type, event.key)
				if event.key.keysym.scancode == .ESCAPE && event.key.state == sdl.PRESSED {
					should_quit = true
				}

				if event.key.keysym.scancode == .RIGHT && event.key.state == sdl.PRESSED {
					pressure_iterations += 1
				}
				if event.key.keysym.scancode == .LEFT  && event.key.state == sdl.PRESSED {
					pressure_iterations -= 1
				}

				pressure_iterations = max(0, pressure_iterations)
			}
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
				for j in 1..<grid_height {
					for i in 1..<grid_width {
						pE := pressure_ping[j][i+1] if i < grid_width-1 else 0.0
						pW := pressure_ping[j][i-1]
						pN := pressure_ping[j+1][i] if j < grid_height-1 else 0.0
						pS := pressure_ping[j-1][i]

						p := (divergence[j][i] + pW + pE + pS + pN) / 4.0
						pressure_pong[j][i] = p
					}
				}
				pressure_ping, pressure_pong = pressure_pong, pressure_ping
			}
		}

		{
			start_timer(&residual_timer)
			defer stop_timer(&residual_timer)

			// Calculate Residual
			for j in 1..<grid_height {
				for i in 1..<grid_width {
					pE := pressure_ping[j][i+1] if i < grid_width-1 else 0.0
					pC := pressure_ping[j][i]
					pW := pressure_ping[j][i-1]
					pN := pressure_ping[j+1][i] if j < grid_height-1 else 0.0
					pS := pressure_ping[j-1][i]

					r := divergence[j][i] + pW + pE + pS + pN - 4 * pC
					residual[j][i] = r
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

			divergence_to_u32 :: proc(div: f32) -> u32 {
				r := u32(255 * clamp(div, 0.0, 1.0))
				g := u32(255 * clamp(-div, 0.0, 1.0))
				b := u32(0)
				a := u32(255)
				return (r << 0) | (g << 8) | (b << 16) | (a << 24) 
			}

			for y in 0..<u32(grid_height) {
				for x in 0..<u32(grid_width) {
					locked_pixels[u32(2*grid_width)*y + x] = velocity_to_u32(velocity_pong[y][x])
				}
				for x in 0..<u32(grid_width) {
					locked_pixels[u32(2*grid_width)*y + x + u32(grid_width)] = divergence_to_u32(divergence[y][x])
				}
			}

			for y in 0..<u32(grid_height) {
				for x in 0..<u32(grid_width) {
					locked_pixels[u32(2*grid_width)*(y + u32(grid_height)) + x] = divergence_to_u32(pressure_ping[y][x])
				}
				for x in 0..<u32(grid_width) {
					locked_pixels[u32(2*grid_width)*(y + u32(grid_height)) + x + u32(grid_width)] = divergence_to_u32(residual[y][x])
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
			render_string(font, renderer, fmt.tprintf("Iterations: %d\x00", pressure_iterations), 0, cursor, text_color); cursor += 16
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

