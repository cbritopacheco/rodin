let SessionLoad = 1
let s:so_save = &g:so | let s:siso_save = &g:siso | setg so=0 siso=0 | setl so=-1 siso=-1
let v:this_session=expand("<sfile>:p")
silent only
silent tabonly
cd ~/Projects/rodin/src/RodinExternal/MMG
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +51 ImplicitDomainMesher2D.h
badd +0 Distancer2D.h
badd +68 ImplicitDomainMesher2D.cpp
badd +801 ~/Projects/rodin/third-party/mmg/src/mmg2d/mmg2d.c
badd +389 ~/Projects/rodin/third-party/mmg/src/mmg2d/API_functions_2d.c
badd +2390 ~/Projects/rodin/third-party/mmg/src/mmg3d/API_functions_3d.c
badd +141 ~/Projects/rodin/third-party/mmg/src/mmg3d/mmg3d.h
badd +64 ~/Projects/rodin/third-party/mmg/src/mmg2d/mmg2d.h
badd +524 ~/Projects/rodin/build/third-party/mmg/include/mmg/mmg2d/libmmgtypes.h
badd +714 ~/Projects/rodin/third-party/mmg/src/common/API_functions.c
badd +113 ~/Projects/rodin/third-party/mmg/src/common/mmgcommon.h
badd +20 ../../../examples/MMG/mmg2d/ImplicitDomainMeshing2D.cpp
badd +809 ~/Projects/rodin/third-party/mmg/src/mmg2d/libmmg2d.c
badd +300 ~/Projects/rodin/build/third-party/mmg/include/mmg/mmg2d/libmmg2d.h
badd +970 ~/Projects/rodin/third-party/mmg/src/mmg2d/mmg2d6.c
badd +312 ~/Projects/rodin/third-party/mmg/src/common/mmg2.c
badd +119 ~/Projects/rodin/third-party/mmg/src/mmg2d/libmmg2d_tools.c
badd +173 ~/Projects/rodin/third-party/mmg/src/common/libmmgcommon.h
badd +184 ~/Projects/rodin/third-party/mmg/src/mmg2d/API_functionsf_2d.c
badd +77 /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include/string.h
argglobal
%argdel
$argadd ImplicitDomainMesher2D.h
set stal=2
tabnew
tabnew
tabrewind
edit ImplicitDomainMesher2D.cpp
argglobal
balt ImplicitDomainMesher2D.h
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let &fdl = &fdl
let s:l = 61 - ((25 * winheight(0) + 17) / 35)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 61
normal! 047|
tabnext
edit ~/Projects/rodin/third-party/mmg/src/mmg2d/libmmg2d_tools.c
argglobal
balt ~/Projects/rodin/third-party/mmg/src/mmg2d/mmg2d.c
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let &fdl = &fdl
let s:l = 111 - ((34 * winheight(0) + 17) / 35)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 111
normal! 05|
tabnext
edit ../../../examples/MMG/mmg2d/ImplicitDomainMeshing2D.cpp
argglobal
balt Distancer2D.h
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let &fdl = &fdl
let s:l = 24 - ((23 * winheight(0) + 17) / 35)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 24
normal! 011|
tabnext 1
set stal=1
if exists('s:wipebuf') && len(win_findbuf(s:wipebuf)) == 0&& getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 shortmess=filnxtToOFAc
let s:sx = expand("<sfile>:p:r")."x.vim"
if filereadable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &g:so = s:so_save | let &g:siso = s:siso_save
set hlsearch
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
