-- Copyright (C) 1999 CAL FP Bank
-- Licensed under Eiffel Forum Freeware License, version 1;
-- (see forum.txt)
--
indexing

	description: "Abstract interface to Fourier Transform"

deferred class FOURIER_OPERATOR

feature

	forward_fft is
		require
			ready: ready_for_forward
		deferred
		ensure
			last_direction = 1
		end

	inverse_fft is
		require
			ready:  ready_for_inverse
		deferred
		ensure
			last_direction = -1
		end

	last_direction: INTEGER
			-- +1 for forward, -1 for inverse, 0 for none


	set_last_direction (llast_direction : like last_direction) is
		require
			corect_values : llast_direction = -1 or llast_direction = 1 or llast_direction = 0
		do
			last_direction := llast_direction
		end


	ready_for_forward, ready_for_inverse: BOOLEAN is
			-- return True, if we are ready to do the transform
		deferred
		end

	is_power_of_two (n: INTEGER): BOOLEAN is
		local
			i: INTEGER
		do
			Result := True
			if i > 1 then
				from i := n
				until (i = 2) or not Result
				loop
					Result := (i \\ 2) = 0
					i := i // 2
				end
			end
		end


end

