# the path to segy file
work_dir = "/Users/wenlei/Desktop/segy";
path_sgy = joinpath(work_dir, "shot_gather.sgy");

# read the text file header
thdr = print_text_header(path_sgy);
