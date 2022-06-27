package ambergin

import sdl "vendor:sdl2"
import "core:math"

Timer :: struct($N: int) {
	counter_frequency: u64,
	previous_counter:  u64,
	
	history:       [N]u64,
	history_index: int,

	sum:         u64,
	squared_sum: u64,

	average: f64,
	std:     f64, // Variation in a single measurement
	ste:     f64, // Variation in the average
}

create_timer :: proc($N: int) -> Timer(N) where N & (N-1) == 0 {
	using timer: Timer(N)

	counter_frequency = sdl.GetPerformanceFrequency()
	previous_counter = sdl.GetPerformanceCounter()

	return timer
}

update_timer :: proc(using timer: ^Timer($N)) {
	current_counter := sdl.GetPerformanceCounter()
	delta_counter := current_counter - previous_counter
	previous_counter = current_counter
	
	old_history_entry := history[history_index & (N-1)]
	new_history_entry := delta_counter
	history[history_index & (N-1)] = new_history_entry

	sum += new_history_entry - old_history_entry
	squared_sum += new_history_entry * new_history_entry - old_history_entry * old_history_entry

	history_index += 1

	num := min(history_index, N)

	average = f64(sum) / f64(num)
	squared_average := f64(squared_sum) / f64(num)
	std = math.sqrt(squared_average - average * average)
	ste = std / math.sqrt(f64(num))

	average *= 1000 / f64(counter_frequency)
	std     *= 1000 / f64(counter_frequency)
	ste     *= 1000 / f64(counter_frequency)
}