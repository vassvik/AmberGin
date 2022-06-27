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

	window_width, window_height: i32 = 1280, 720
	window := sdl.CreateWindow("Test SDL", 50, 50, window_width, window_height, {.SHOWN})
	defer sdl.DestroyWindow(window)
	assert(window != nil)

	renderer := sdl.CreateRenderer(window, -1, {.ACCELERATED})
	//renderer := sdl.CreateRenderer(window, -1, {.ACCELERATED, .PRESENTVSYNC})
	defer sdl.DestroyRenderer(renderer)
	assert(renderer != nil)

	tex := sdl.CreateTexture(renderer, u32(sdl.PixelFormatEnum.ABGR8888), .STREAMING, window_width, window_height)
	defer sdl.DestroyTexture(tex)
	assert(tex != nil)

	frame_timer := create_timer(128)

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


	velocity_ping := make_2D([2]f32, window_width, window_height)
	velocity_pong := make_2D([2]f32, window_width, window_height)

	divergence := make_2D(f32, window_width, window_height)
	pressure_ping := make_2D(f32, window_width, window_height)
	pressure_pong := make_2D(f32, window_width, window_height)

	should_quit := false
	for i := u32(0); !should_quit; i += 1 {
		update_timer(&frame_timer)

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

		// Stream texture
		@static time_copy := 0.0;
		{
			// Get texture memory
			locked_pixels_raw: rawptr
			pitch: i32
			sdl.LockTexture(tex, nil, &locked_pixels_raw, &pitch)
			defer sdl.UnlockTexture(tex)
			assert(locked_pixels_raw != nil)
			assert(pitch == window_width * size_of([4]u8))

			// Write to pixels
			locked_pixels := mem.slice_ptr(cast(^u32)locked_pixels_raw, int(window_width) * int(window_height))

			u8x4_to_u32 :: proc(r, g, b, a: u8) -> u32 {
				return (u32(r) << 0) | (u32(g) << 8) | (u32(b) << 16) | (u32(a) << 24) 
			}

			t1 := sdl.GetPerformanceCounter()
			for y in 0..<u32(window_height) {
				for x in 0..<u32(window_width) {
					locked_pixels[u32(window_width)*y + x] = u8x4_to_u32(u8(x + 2*i), u8(y + i), 0, 255)
				}
			}
			t2 := sdl.GetPerformanceCounter()
			dt := f64(t2 - t1) / f64(sdl.GetPerformanceFrequency())

			omega := 0.02
			time_copy = time_copy * (1.0 - omega) + dt * omega
		}

		sdl.RenderClear(renderer);
		sdl.RenderCopy(renderer, tex, nil, nil);

		// Render Text		
		text_color := sdl.Color{255, 255, 255, 255};
		render_string(font, renderer, fmt.tprintf("Frame: %d", i), 0, 0, text_color)
		render_string(font, renderer, fmt.tprintf("Timer: %f +/- %f (%f) ms", frame_timer.average, frame_timer.std, frame_timer.ste), 0, 16, text_color)
		render_string(font, renderer, fmt.tprintf("Array Write: %f ms", 1000*time_copy), 0, 32, text_color)
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
