using BinDeps

@BinDeps.setup

# set to true to support intel fortran compiler
useIntelFortran = false

# construct absolute path to the dependency directory
depsdir  = splitdir(Base.source_path())[1]
builddir = joinpath(depsdir,"builds")
srcdir   = joinpath(depsdir,"src")

# print the directory
println("=== Building spmatvec ===")
println("depsdir  = $depsdir")
println("builddir = $builddir")
println("srcdir   = $srcdir")
println("useIntel = $useIntelFortran")

# create directory to save the binary
if !isdir(builddir)
	 println("creating build directory")
	 mkdir(builddir)
	 if !isdir(builddir)
		  error("Could not create build directory")
	 end
end

# compile source code for linux machine
@static if Sys.islinux()

	path_src = joinpath(srcdir,"spmatvec.f90")
	path_bin = joinpath(builddir,"spmatvec")

	@build_steps begin
		if useIntelFortran
			 run(`ifort -O3 -xHost -fPIC -fpp -openmp -integer-size 64 -diag-disable=7841 -shared  $path_src -o $path_bin`)
		else
			 println("fortran version")
			 run(`gfortran --version`)
			 run(`gfortran -O3 -fPIC -cpp -fopenmp -fdefault-integer-8 -shared  $path_src -o $path_bin`)
		end
	end
end

# compile source code for Mac
@static if Sys.isapple()

  path_src = joinpath(srcdir,"spmatvec.f90")
	path_bin = joinpath(builddir,"spmatvec")

	@build_steps begin
		if useIntelFortran
			 run(`ifort -O3 -xHost -fPIC -fpp -openmp -integer-size 64 -diag-disable=7841 -shared  $path_src -o $path_bin`)
		else
			 println("fortran version")
			 run(`gfortran --version`)
			 run(`gfortran -O3 -fPIC -cpp -fopenmp -fdefault-integer-8 -shared  $path_src -o $path_bin`)
		end
	end
end
