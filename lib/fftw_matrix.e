-- Copyright (C) 1999 CAL FP Bank
-- Licensed under Eiffel Forum Freeware License, version 1;
-- (see forum.txt)
--
indexing

	description: "Matrix for interfacing with FFTW (http://www.fftw.org)%
                 % Matrices are stored in row-major order."
	author: "Richie Bielak"
	date: "1999/12/07"



class FFTW_MATRIX

inherit

	EXCEPTIONS
		redefine
			copy, is_equal
		end

creation

	set_spec_from_vector,
	set_spec_zero,
	make

feature -- creation

	set_spec_from_vector (v: ARRAY [DOUBLE]; big_row, big_col, little_row, little_col: INTEGER) is
		require
			valid_vector: v /= Void
			consitent_row: (big_row >= little_row) and (little_row > 0)
			consistent_column: (big_col >= little_col) and (little_col > 0)
			consistent_with_vector: (little_col * little_row) = v.count
		local
			r, c, i: INTEGER
			to_area: TO_SPECIAL [DOUBLE]
		do
			make (big_row, big_col)
			-- load the vector into the matrix
			i := 1
			from r := 1 until r > little_row
			loop
				from c := 1 until c > little_col
				loop
					-- put (v @ i, r, c)
					area.put (v @ i, (r - 1) * col_size + c - 1)
					c := c + 1
					i := i + 1
				end
				r := r + 1
			end
		end

	set_spec_zero (row, col: INTEGER) is
			-- create a matrix with all zero entries
		require
			row_valid: row > 0
			col_valid: col > 0
		local
			i: INTEGER
		do
			make (row, col)
			from i := 0 until i >= row_size * col_size
			loop
				area.put (0.0, i)
				i := i + 1
			end
		end

	make (row, col: INTEGER) is
			-- create a matrix and take default initialization
		require
			row_valid: row > 0
			col_valid: col > 0
		local
			to_area: TO_SPECIAL [DOUBLE]
		do
			row_size := row
			col_size := col
			-- allocate an area
			-- !!area.make (row_size * col_size)
			!!to_area.make_area (row_size * col_size)
			area := to_area.area
		end

feature -- size and access

	row_size: INTEGER

	col_size: INTEGER

	item (row, col: INTEGER): DOUBLE is
		require
			row_valid: (row > 0) and (row <= row_size)
			col_valid: (col > 0) and (col <= col_size)
		do
			Result := area.item ((row - 1) * col_size + col - 1)
		end

	put (new_item: DOUBLE; row, col: INTEGER) is
		require
			row_valid: (row > 0) and (row <= row_size)
			col_valid: (col > 0) and (col <= col_size)
		do
			area.put (new_item, (row - 1) * col_size + col - 1)
		ensure
			added: new_item = item (row, col)
		end

feature -- conversions

	to_vector (row, col: INTEGER): ARRAY [DOUBLE] is
		require
			row_valid: (row > 0) and (row <= row_size)
			col_valid: (col > 0) and (col <= col_size)
		do
			!!Result.make (1, row*col)
			fill_vector (Result, row, col)
		ensure
			result_size: Result.count = row * col	
		end

	fill_vector (v: ARRAY [DOUBLE]; row, col: INTEGER) is
		require
			valid_vector: (v /= Void) and then (v.count = row * col)
			row_valid: (row > 0) and (row <= row_size)
			col_valid: (col > 0) and (col <= col_size)
		local
			i, r, c: INTEGER
			it: DOUBLE
		do
			i := 1
			from r := 1 until r > row loop
				from c := 1 until c > col loop
					-- access area directly to speed things up
					it := area.item ((r - 1) * col_size + c - 1)
					v.put (it, i)
					i := i + 1
					c := c + 1
				end
				r := r + 1
			end
		end

	fill_and_scale_vector (v: ARRAY [DOUBLE]; factor: DOUBLE; row, col: INTEGER) is
		require
			valid_vector: (v /= Void) and then (v.count = row * col)
			row_valid: (row > 0) and (row <= row_size)
			col_valid: (col > 0) and (col <= col_size)
		local
			i, r, c: INTEGER
			it: DOUBLE
		do
			i := 1
			from r := 1 until r > row loop
				from c := 1 until c > col loop
					-- access area directly to speed things up
					it := area.item ((r - 1) * col_size + c - 1)
					v.put (it * factor, i)
					i := i + 1
					c := c + 1
				end
				r := r + 1
			end
		end


feature -- operations

	real_element_multiply_by (other: like Current) is
			-- multiply elementwise
		require
			other.row_size = row_size
			other.col_size = col_size
		local
			i: INTEGER
			temp: DOUBLE
		do
			from i := 0 until i >= area.count
			loop
				temp := other.area.item (i) * area.item (i)
				area.put (temp, i)
				i := i + 1
			end
		ensure
			-- result placed in current, other unchanged
		end


	complex_element_multiply_by (other: like Current) is
			-- Treat pairs of entries in the matrix rows as complex 
			-- numbers and multiply the complex elements of the two 
			-- matrices placing the result in Current
		require
			compatible_size: (other.row_size = row_size) and (other.col_size = col_size)
			valid_columns: (col_size \\ 2) = 0
		local
			x, y, a1, a2, b1, b2: DOUBLE
			r, c: INTEGER
		do
			-- this code should be re-written in C, if it turns out 
			-- to be used a lot.
			from r := 1 until r > row_size loop
				from c := 1 until c > col_size loop
					-- Use area and save temporary results to improve 
					-- the speed
					a1 := item (r,c)   -- re (item)
					a2 := item (r,c+1) -- im (item)
					b1 := other.item (r,c) -- re (other.item)
					b2 := other.item (r,c+1) -- im (other.item)
					-- new real part
					-- x := item (r,c)*other.item (r,c) - item(r,c+1)*other.item(r,c+1)
					x := a1 * b1 - a2 * b2
					-- new imaginary part
					-- y := item (r,c+1)*other.item (r,c) + item (r,c)*other.item(r,c+1)
					y := a2 * b1 + a1 * b2
					-- save in current matrix
					-- put (x, r, c)
					area.put (x, (r - 1) * col_size + c - 1)
					-- put (y, r, c + 1)
					area.put (y, (r - 1) * col_size + c)
					c := c + 2
				end
				r := r + 1
			end
		end


	element_add (other: like Current) is
			-- add elements of other to current
		require
			other.row_size = row_size
			other.col_size = col_size
		local
			i: INTEGER
			temp: DOUBLE
		do
			from i := 0 until i >= area.count
			loop
				temp := other.area.item (i) + area.item (i)
				area.put (temp, i)
				i := i + 1
			end
		end

	scalar_multiply (scalar: DOUBLE) is
			-- muliply all item by a scalar
		local
			i: INTEGER
		do
			from i := 0 until i >= area.count
			loop
				area.put (area.item (i) * scalar, i)
				i := i + 1
			end
		end

	conjugate is
			-- replace the matrix by its conjugate (real element 
			-- remains equal, imaginary are * -1)
		require
			col_size_even: (col_size \\ 2) = 0
		local
			x: DOUBLE
			r, c: INTEGER
		do
			from r := 1 until r > row_size loop
				from c := 2 until c > col_size loop
					x := item (r, c)
					put (-x, r, c)
					c := c + 2
				end
				r := r + 1
			end
		end		
		
	dump is
		local
			i, j : INTEGER
		do
			print ("FFTW Matrix: ")
			print (row_size);
			print (" x ")
			print (col_size)
			print ("%N")
			from i := 1
			until (i > row_size)
			loop
				print ("[ ")
				from j := 1
				until j > col_size
				loop
					print (item (i, j))
					print (" ")
					j := j + 1
				end
				print (" ]%N")
				i := i + 1
			end
		end


feature -- standard features

	copy (other: like Current) is
		do
			row_size := other.row_size
			col_size := other.col_size
			area := clone (other.area)
		end

	is_equal (other: like Current): BOOLEAN is
		do
			Result := (row_size = other.row_size) and (col_size = other.col_size)
			if Result then
				Result := area.is_equal (other.area)
			end
		end

	quasi_is_equal (other: like Current): BOOLEAN is
		require
			other_not_void: other /= Void
		local
			i: INTEGER
			abs: DOUBLE
		do
			Result := (row_size = other.row_size) and (col_size = other.col_size)
			if Result then
				from i := 0 until (i >= area.count) or not Result
				loop
					abs := (area.item (i) - other.area.item (i))
					if abs < 0.0 then 
						abs := - abs
					end
					Result := abs <= 1.0e-15
					i := i + 1
				end
			end
		end

	to_c: ANY is
		do
			Result := area
		end

feature {FFTW_MATRIX} -- implementation

	area: SPECIAL [DOUBLE]

end
