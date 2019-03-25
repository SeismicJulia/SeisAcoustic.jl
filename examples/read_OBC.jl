using SeisAcoustic

# the path to segy file
work_dir = "/home/wgao1/australia_OBN";
path_sgy = joinpath(work_dir, "hydro.segy");

# read text file header
thdr = read_text_header(path_sgy);

# read file header
fhdr = read_file_header(path_sgy);

# read all the trace header in one segy file
traces_header = read_all_traces_header(path_sgy);

# number of traces
tmp = filesize(path_sgy);
num_traces = convert(Int64, (tmp - 3600) / (240+sizeof(Float32)*3000))

# read all the traces header
h = init_TraceHeader()
