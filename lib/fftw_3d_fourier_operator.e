-- Copyright (C) 1999 CAL FP Bank
-- Licensed under Eiffel Forum Freeware License, version 1;
-- (see forum.txt)
--
indexing

	description: "3D real fourier transform operator. Uses the FFTW %
	              %libraries. See: http://www.fftw.org"
	author: "Richie Bielak"
	date: "December 14th, 1999"

class FFTW_3D_FOURIER_OPERATOR

inherit

	FOURIER_OPERATOR

	FFTW_EXTERNALS

	MEMORY
		redefine
			dispose
		end

creation

	make_both,
	make_forward_plan_only,
	make_inverse_plan_only

feature -- creation

	make_both (lrows, lcols, ldepth: INTEGER) is
		require
			valid_rows: lrows > 0
			valid_cols: lcols > 0
			valid_depth: ldepth > 0
		do
			make_forward_plan_only (lrows, lcols, ldepth)
			make_inverse_plan_only (lrows, lcols, ldepth)
		ensure
			size_set: (rows = lrows) and (cols = lcols) and (depth = ldepth)
			iplan_valid: inverse_plan /= default_pointer
			fplan_valid: forward_plan /= default_pointer
		end

	make_forward_plan_only (lrows, lcols, ldepth: INTEGER) is
		require
			valid_rows: lrows > 0
			valid_cols: lcols > 0
			valid_depth: ldepth > 0
		do
			rows := lrows
			cols := lcols
			depth := ldepth
			if forward_plan /= default_pointer then
				rfftwnd_destroy_plan (forward_plan)
			end
			forward_plan := rfftw3d_create_plan (rows, cols, depth, -1, 0)
		ensure
			size_set: (rows = lrows) and (cols = lcols) and (depth = ldepth)
			fplan_valid: forward_plan /= default_pointer
		end

	make_inverse_plan_only (lrows, lcols, ldepth: INTEGER) is
		require
			valid_rows: lrows > 0
			valid_cols: lcols > 0
			valid_depth: ldepth > 0
		do
			rows := lrows
			cols := lcols
			depth := ldepth
			if inverse_plan /= default_pointer then
				rfftwnd_destroy_plan (inverse_plan)
			end
			inverse_plan := rfftw3d_create_plan (rows, cols, depth, 1, 0)
		ensure
			size_set: (rows = lrows) and (cols = lcols) and (depth = ldepth)
			iplan_valid: inverse_plan /= default_pointer
		end

feature -- queries

	rows: INTEGER

	cols: INTEGER

	depth: INTEGER

	ready_for_forward: BOOLEAN is
		do
			Result := (forward_plan /= default_pointer) 
				and then real_tensor /= Void
				and then complex_tensor /= Void
		end

	ready_for_inverse: BOOLEAN is
		do
			Result := inverse_plan /= default_pointer
				and then real_tensor /= Void
				and then complex_tensor /= Void
		end


feature -- operations

	real_tensor: FFTW_3DTENSOR

	set_real (t: like real_tensor) is
		require
			t_not_void: t /= Void
			consistent_size: (t.row_size = rows) and (t.col_size = cols) and (t.depth_size = depth)
		do
			real_tensor := t
		ensure
			set: real_tensor = t
		end

	complex_tensor: like real_tensor

	set_complex (t: like real_tensor) is
		require
			t_not_void: t /= Void
			consistent_size: (t.row_size = rows) and (t.col_size = cols)
			consistent_depth: t.depth_size  = 2 *((depth // 2) + 1)
		do
			complex_tensor := t
		ensure
			set: complex_tensor = t
		end

	forward_fft is
		do
			rfftwnd_one_real_to_complex (forward_plan, $(real_tensor.to_c), $(complex_tensor.to_c))
			last_direction := 1
		end

	inverse_fft is
		do
			rfftwnd_one_complex_to_real (inverse_plan, $(complex_tensor.to_c), $(real_tensor.to_c))
			last_direction := -1
		end

feature {NONE} -- implementation

	forward_plan: POINTER

	inverse_plan: POINTER

	dispose is
		do
			if forward_plan /= default_pointer then
				rfftwnd_destroy_plan (forward_plan)
			end
			if inverse_plan /= default_pointer then
				rfftwnd_destroy_plan (inverse_plan)
			end
		end

invariant

	plan_defined: (forward_plan /= default_pointer) or (inverse_plan /= default_pointer)
	consistent_dimension: (rows > 0) and (cols > 0) and (depth > 0)

end
