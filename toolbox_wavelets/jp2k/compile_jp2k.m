cd src

% remove switch -Dcompil_vcc is you are under linux/unix

mex -g -O -Dcompil_vcc -DWIN32 -outdir ../ -output jp2_class                     jp2_codec.c bio.c dwt.c j2k.c mct.c pi.c t2.c tgt.c cio.c fix.c int.c mqc.c t1.c tcd.c image_jp2.c liw_jp2_dll.c liw_error.c

mex -g -O -Dcompil_vcc -DWIN32 -DENCODE_ONLY -outdir ../ -output jp2_codec_encod jp2_codec.c bio.c dwt.c j2k.c mct.c pi.c t2.c tgt.c cio.c fix.c int.c mqc.c t1.c tcd.c image_jp2.c liw_jp2_dll.c liw_error.c

cd ..