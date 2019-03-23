using LinearAlgebra, SeisAcoustic

# the path to segy file
work_dir = "/Users/wenlei/Desktop/segy";
path_sgy = joinpath(work_dir, "shot_gather.sgy");

# # read the binary
# fid = open(path_sgy, "r");
# thdr = read(fid, 3200);
# close(fid);
#
# #convert to ascii
# thdr1 = SeisAcoustic.conversion_between_edbcb_ascii(thdr, flag="edbcb2ascii")
# thdr2 = convert(Vector{Char}, thdr1)
#
# # print the result
# for i = 1 : 40
#     il = (i-1)*80+1
#     iu = il+80-1
#     println(join(thdr2[il:iu]))
# end
#
# # return back to edbcb
# thdr3 = convert(Vector{UInt8}, thdr2);
# thdr4 = SeisAcoustic.conversion_between_edbcb_ascii(thdr3, flag="ascii2edbcb")
# norm(thdr4 - thdr)
#
# read the text file header
thdr = print_text_header(path_sgy);
#
# # test write a segy file
path_tmp = joinpath(work_dir, "tmp.sgy");
put_text_header(path_tmp, thdr);
#
# # read the text header out of generated file
# thdr1 = print_text_header(path_tmp);

# test read the file header
fhdr = read_file_header(path_sgy);

# test write the file header
put_file_header(path_tmp, fhdr);


fhdr1 = read_file_header(path_tmp);
b = diff_file_header(fhdr, fhdr1);

fhdr2 = FileHeader(fhdr);
# test they are equal to each other
function diff_file_header(fhdr1, fhdr2)

    r = 0.0
    for field in fieldnames(FileHeader)
        r = r + abs2(getfield(fhdr1, field) - getfield(fhdr2, field))
    end
    return r
end
