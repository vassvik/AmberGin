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

	window := sdl.CreateWindow("Test SDL", 50, 50, 1280, 720, {.SHOWN})
	defer sdl.DestroyWindow(window)
	assert(window != nil)

	renderer := sdl.CreateRenderer(window, -1, {.ACCELERATED})
	//renderer := sdl.CreateRenderer(window, -1, {.ACCELERATED, .PRESENTVSYNC})
	defer sdl.DestroyRenderer(renderer)
	assert(renderer != nil)

	w, h: i32 = 256, 256
	tex := sdl.CreateTexture(renderer, u32(sdl.PixelFormatEnum.ABGR8888), .STREAMING, w, h)
	defer sdl.DestroyTexture(tex)
	assert(tex != nil)


	frame_timer := create_timer(128)

	previous_qpc: u64
	qpf := sdl.GetPerformanceFrequency()
	
	qpc_history: [32]u64
	qpc_average: u64
	pqc_counter := 0;

	should_quit := false
	for i := 0; !should_quit; i += 1 {
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
		{
			// Get texture memory
			locked_pixels_raw: rawptr
			pitch: i32
			sdl.LockTexture(tex, nil, &locked_pixels_raw, &pitch)
			defer sdl.UnlockTexture(tex)
			assert(locked_pixels_raw != nil)
			assert(pitch == w * size_of([4]u8))

			// Write to pixels
			locked_pixels := mem.slice_ptr(cast(^[4]u8)locked_pixels_raw, int(w * h))
			for y in 0..<h {
				for x in 0..<w {
					locked_pixels[w*y + x] = {u8(int(x) + 2*i), u8(int(y) + i), 0, 255}
				}
			}
		}

		sdl.RenderClear(renderer);
		sdl.RenderCopy(renderer, tex, nil, nil);

		// Render Text		
		text_color := sdl.Color{255, 255, 255, 255};
		render_string(font, renderer, fmt.tprintf("Frame: %d", i), 0, 0, text_color)
		render_string(font, renderer, fmt.tprintf("Timer: %f +/- %f (%f)", frame_timer.average, frame_timer.std, frame_timer.ste), 0, 32, text_color)
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
