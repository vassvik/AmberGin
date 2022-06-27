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

	window_width, window_height: i32 = 2*512, 512
	window := sdl.CreateWindow("Test SDL", 50, 50, window_width, window_height, {.SHOWN})
	defer sdl.DestroyWindow(window)
	assert(window != nil)

	//renderer := sdl.CreateRenderer(window, -1, {.ACCELERATED})
	renderer := sdl.CreateRenderer(window, -1, {.ACCELERATED, .PRESENTVSYNC})
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

	grid_width, grid_height: i32 = 64, 64
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
	pressure_ping := make_2D(f32, grid_width, grid_height)
	pressure_pong := make_2D(f32, grid_width, grid_height)

	display_texture := sdl.CreateTexture(renderer, u32(sdl.PixelFormatEnum.ABGR8888), .STREAMING, 2*grid_width, grid_height)
	defer sdl.DestroyTexture(display_texture)
	assert(display_texture != nil)


	frame_timer := create_timer(128)
	start_timer(&frame_timer)

	display_write_timer := create_timer(128)
	divergence_timer := create_timer(128)

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
			}
		}

		{
			start_timer(&divergence_timer)
			defer stop_timer(&divergence_timer)
			// Calculate divergence
			for j in 0..<grid_height {
				for i in 0..<grid_width {
					if i == 0 || j == 0 do continue

					lower_left := velocity_ping[j-1][i-1]
					lower_right := velocity_ping[j-1][i]
					upper_left := velocity_ping[j][i-1]
					upper_right := velocity_ping[j][i]

					vE := 0.5 * (lower_right.x + upper_right.x)
					vW := 0.5 * (lower_left.x + upper_left.x)
					vN := 0.5 * (upper_left.y + upper_right.y)
					vS := 0.5 * (lower_left.y + lower_right.y)

					divergence[j][i] = (vE - vW) + (vN - vS)
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
			locked_pixels := mem.slice_ptr(cast(^u32)locked_pixels_raw, int(2*grid_width) * int(grid_height))

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
					locked_pixels[u32(2*grid_width)*y + x] = velocity_to_u32(velocity_ping[y][x])
				}
				for x in 0..<u32(grid_width) {
					locked_pixels[u32(2*grid_width)*y + x + u32(grid_width)] = divergence_to_u32(divergence[y][x])
				}
			}
		}

		sdl.RenderClear(renderer);
		sdl.RenderCopy(renderer, display_texture, nil, nil);

		// Render Text		
		cursor := i32(0)
		text_color := sdl.Color{255, 255, 255, 255};
		render_string(font, renderer, fmt.tprintf("Frame: %d", frame), 0, cursor, text_color); cursor += 16
		render_string(font, renderer, fmt.tprintf("Frame Timer: %f +/- %f (%f) ms", frame_timer.average, frame_timer.std, frame_timer.ste), 0, cursor, text_color); cursor += 16
		render_string(font, renderer, fmt.tprintf("Display Write Timer: %f +/- %f (%f) ms", display_write_timer.average, display_write_timer.std, display_write_timer.ste), 0, cursor, text_color); cursor += 16
		render_string(font, renderer, fmt.tprintf("Divergence Timer: %f +/- %f (%f) ms", divergence_timer.average, divergence_timer.std, divergence_timer.ste), 0, cursor, text_color); cursor += 16
		sdl.RenderPresent(renderer);
	}
}


render_string :: proc(font: ^sdl_ttf.Font, renderer: ^sdl.Renderer, str: string, x, y: i32, color: sdl.Color) {
	cstr := strings.unsafe_string_to_cstring(str)

	text := sdl_ttf.RenderText_Blended(font, cstr, color)
	assert(text != nil)

	text_tex := sdl.CreateTextureFromSurface(renderer, text)
	assert(text_tex != nil)

	rect := sdl.Rect{x, y, text.w, text.h}
	sdl.RenderCopy(renderer, text_tex, nil, &rect);
}
