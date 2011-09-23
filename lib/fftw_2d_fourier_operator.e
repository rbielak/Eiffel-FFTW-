-- Copyright (C) 1999 CAL FP Bank
-- Licensed under Eiffel Forum Freeware License, version 1;
-- (see forum.txt)
--
indexing

	description: "2D real fourier transform operator. Uses the FFTW %
                 %libraries. See: http://www.fftw.org"
	author: "Richie Bielak"
	date: "December 10th, 1999"

class FFTW_2D_FOURIER_OPERATOR

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

	make_both (lrows, lcols: INTEGER) is
			-- create plans for both forward and inverse transformations
		require
			valid_rows: lrows > 0
			valid_cols: lcols > 0
		do
			make_forward_plan_only (lrows, lcols)
			make_inverse_plan_only (lrows, lcols)
		ensure
			size_set: (rows = lrows) and (cols = lcols)
			iplan_valid: inverse_plan /= default_pointer
			fplan_valid: forward_plan /= default_pointer
		end

	make_forward_plan_only (lrows, lcols: INTEGER) is
			-- create plans for forward transformation only
		require
			valid_rows: lrows > 0
			valid_cols: lcols > 0
		do
			rows := lrows
			cols := lcols
			if forward_plan /= default_pointer then
				rfftwnd_destroy_plan (forward_plan)
			end
			forward_plan := rfftw2d_create_plan (rows, cols, -1, 0)
		ensure
			size_set: (rows = lrows) and (cols = lcols)
			plan_valid: forward_plan /= default_pointer
		end

	make_inverse_plan_only (lrows, lcols: INTEGER) is
			-- create plans for forward transformation only
		require
			valid_rows: lrows > 0
			valid_cols: lcols > 0
		do
			rows := lrows
			cols := lcols
			if inverse_plan /= default_pointer then
				rfftwnd_destroy_plan (inverse_plan)
			end
			inverse_plan := rfftw2d_create_plan (rows, cols, 1, 0) 
		ensure
			size_set: (rows = lrows) and (cols = lcols)
			plan_valid: inverse_plan /= default_pointer
		end

feature -- queries

	rows: INTEGER
			-- number of rows in the input matrix

	cols: INTEGER
			-- number of cols in the input matrix


	ready_for_forward: BOOLEAN is
			-- Forward transform takes a real "rows x cols" matrix and 
			-- transforms it into a complex matrix with rows x 
			-- (cols//2 + 1) complex numbers. This means that the 
			-- actual column index runs from 1 to 2 * ((cols//2) + 1) 
			-- since each complex number has Re and Im part.
		do
			Result := (forward_plan /= default_pointer) 
				and then real_matrix /= Void
				and then complex_matrix /= Void
		end

	ready_for_inverse: BOOLEAN is
			-- Inverse transform takes a complex matrix of size rows x ((cols // 2) + 1)
			-- which means that the column index actually runs from 1 
			-- to 2 * ((cols // 2) + 1). The output matrix is real 
			-- and sized rows x cols.
		do
			Result := (inverse_plan /= default_pointer) 
				and then complex_matrix /= Void
				and then real_matrix /= Void
		end

feature -- operations

	real_matrix: FFTW_MATRIX

	set_real (m: like real_matrix) is
		require
			not_void: m /= Void
			correct_size: (m.row_size = rows) and (m.col_size = cols)
		do
			real_matrix := m
		ensure
			set: real_matrix = m
		end


	complex_matrix: FFTW_MATRIX 

	set_complex (m: like real_matrix) is
			-- The complex matrix, which is the output of forward and 
			-- input to the inverse transform has size:
			-- rows x ((cols // 2) +  1).
			-- Which means that the column index actually runs from 1 
			-- to 2 * ((cols // 2) + 1) in order to hold the real and 
			-- imaginary parts of the numbers

		require
			not_void: m /= Void
			correct_size: (m.row_size = rows) and (m.col_size = 2 * (cols//2 + 1))
		do
			complex_matrix := m
		ensure
			set: complex_matrix = m
		end


	forward_fft is
		do
			rfftwnd_one_real_to_complex (forward_plan, $(real_matrix.to_c), $(complex_matrix.to_c))
			last_direction := 1
		end

	inverse_fft is
		do
			rfftwnd_one_complex_to_real (inverse_plan, $(complex_matrix.to_c), $(real_matrix.to_c))
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
	consistent_dimension: (rows > 0) and (cols > 0) 

end
