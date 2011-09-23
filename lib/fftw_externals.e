class FFTW_EXTERNALS

feature {NONE}

	rfftw_create_plan (n: INTEGER; direction: INTEGER; flags: INTEGER): POINTER is
			-- direction: -1 -> forward, +1 inverse
		external "C"
		end

	rfftw2d_create_plan (r, c: INTEGER; dir: INTEGER; flags: INTEGER): POINTER is
		external "C"
		end

	rfftw3d_create_plan (x, y, z: INTEGER; dir: INTEGER flags: INTEGER): POINTER is
		external "C"
		end

	rfftw_destroy_plan (pp: POINTER) is
		external "C"
		end

	rfftwnd_destroy_plan (p: POINTER) is
		external "C"
		end

	rfftw_one (plan: POINTER; in_v: POINTER; out_v: POINTER) is
		external "C"
		end
	
	rfftwnd_one_real_to_complex (plan: POINTER; in_m: POINTER; out_m: POINTER) is
		external "C"
		end

	rfftwnd_one_complex_to_real (plan: POINTER; in_m: POINTER; out_m: POINTER) is
		external "C"
		end


end
