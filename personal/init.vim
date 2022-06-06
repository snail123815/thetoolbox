" Basics {{{
filetype plugin indent on
syntax on
"}}}

" Settings {{{
:set number
:augroup numbertoggle
:  autocmd!
:  autocmd BufEnter,FocusGained,InsertLeave,WinEnter * if &nu && mode() != "i" | set rnu   | endif
:  autocmd BufLeave,FocusLost,InsertEnter,WinLeave   * if &nu                  | set nornu | endif
:augroup END
set backspace=indent,eol,start    " Backspace everything in insert mode
set history=50
set ruler
set showcmd
set incsearch
set hlsearch
set tabstop=4
set shiftwidth=4
set softtabstop=4
set autoindent
set formatoptions=tcroql          " Auto-wrap comments
set wildmenu                      " Display matches in command-line mode
set expandtab                     " Prefer spaces over tabs in general
set hidden                        " Prefer hiding over unloading buffers
set path=.,**                     " Relative to current file and everything under :pwd
setl wildignore=**/node_modules/**,**/dist/**,*.pyc
set tags=./tags;,tags;            " Find tags relative to current file and directory
set t_BE=                         " Disable bracketed paste mode
set nolangremap
set backup
set backupdir-=.
set backupdir^=~/.nvimbk

" IMPORTANT: :help Ncm2PopupOpen for more information
set completeopt=noinsert,menuone,noselect
"}}}


" Status line {{{
function! GitBranch()
  return system("git rev-parse --abbrev-ref HEAD 2>/dev/null | tr -d '\n'")
endfunction

function! StatuslineGit()
  let l:branchname = GitBranch()
  return strlen(l:branchname) > 0?'  '.l:branchname.' ':''
endfunction

set laststatus=2
set statusline=
set statusline+=%n\ \|
set statusline+=%{StatuslineGit()}
set statusline+=\ %f%m%r%h%w
set statusline+=%=
set statusline+=\ %y
set statusline+=\ %p%%
set statusline+=\ %l:%c
set statusline+=\ 
"}}}

" When git is using nvim to commit (I guess) {{{
autocmd BufRead * autocmd FileType <buffer> ++once
      \ if &ft !~# 'commit\|rebase' && line("'\"") > 1 && line("'\"") <= line("$") | exe 'normal! g`"' | endif
"}}}

" Vim-plug {{{
let data_dir = has('nvim') ? stdpath('data') . '/site' : '~/.vim'
if empty(glob(data_dir . '/autoload/plug.vim'))
  silent execute '!curl -fLo '.data_dir.'/autoload/plug.vim --create-dirs  https://raw.githubusercontent.com/junegunn/vim-plug/master/plug.vim'
  autocmd VimEnter * PlugInstall --sync | source $MYVIMRC
endif

call plug#begin(has('nvim') ? stdpath('data') . '/plugged' : '~/.local/neovim/plugged')

" Declare the list of plugins.

Plug 'ncm2/ncm2'
Plug 'ncm2/ncm2-jedi'
Plug 'ncm2/ncm2-bufword'
Plug 'ncm2/ncm2-path'

Plug 'ojroques/vim-oscyank'
Plug 'tpope/vim-sensible'
Plug 'tpope/vim-surround'
Plug 'roxma/nvim-yarp'
Plug 'preservim/tagbar'
Plug 'junegunn/seoul256.vim'
Plug 'godlygeek/tabular'
Plug 'plasticboy/vim-markdown'
Plug 'mzlogin/vim-markdown-toc'

" List ends here. Plugins become visible to Vim after this call.
call plug#end()

" Plugin settings {

" 'ncm2/ncm2'
" enable ncm2 for all buffers
autocmd BufEnter * call ncm2#enable_for_buffer()
" wrap existing omnifunc
" Note that omnifunc does not run in background and may probably block the
" editor. If you don't want to be blocked by omnifunc too often, you could
" add 180ms delay before the omni wrapper:
"  'on_complete': ['ncm2#on_complete#delay', 180,
"               \ 'ncm2#on_complete#omni', 'csscomplete#CompleteCSS'],
au User Ncm2Plugin call ncm2#register_source({
    \ 'name' : 'css',
    \ 'priority': 9,
    \ 'subscope_enable': 1,
    \ 'scope': ['css','scss'],
    \ 'mark': 'css',
    \ 'word_pattern': '[\w\-]+',
    \ 'complete_pattern': ':\s*',
    \ 'on_complete': ['ncm2#on_complete#omni', 'csscomplete#CompleteCSS'],
    \ })
" IMPORTANT: python3 path
let g:python3_host_prog = '$HOME/apps/micromamba/bin/python3'

" 'plasticboy/vim-markdown'
let g:vim_markdown_toc_autofit = 1
let g:vim_markdown_math = 1

" 'mzlogin/vim-markdown-toc'
let g:vmt_dont_insert_fence = 1
let g:vmt_list_item_char = "-"

" 'ojroques/vim-oscyank'
" Yank to system clipboard
augroup Yank
  autocmd!
  autocmd TextYankPost * if v:event.operator is 'y' && v:event.regname is '' | execute 'OSCYankReg "' | endif
augroup END

"}
"}}}


" Key mappings {{{
" CTRL-C doesn't trigger the InsertLeave autocmd . map to <ESC> instead.
inoremap <c-c> <ESC>

" When the <Enter> key is pressed while the popup menu is visible, it only
" hides the menu. Use this mapping to close the menu and also start a new
" line.
inoremap <expr> <CR> (pumvisible() ? "\<c-y>\<cr>" : "\<CR>")

"}}}


" FileType specific {{{
augroup filetype_md
	autocmd!
	:autocmd FileType markdown setlocal tabstop=2 shiftwidth=2 softtabstop=2 expandtab autoindent
	:autocmd FileType markdown noremap <silent> k gk
	:autocmd FileType markdown noremap <silent> j gj
	:autocmd FileType markdown noremap <silent> 0 g0
	:autocmd FileType markdown noremap <silent> $ g$
augroup end

augroup filetype_py
    autocmd!
    :autocmd FileType python setlocal shiftwidth=4 tabstop=4 softtabstop=4 noexpandtab autoindent smartindent
    :autocmd FileType python setlocal colorcolumn=80
    :autocmd FileType python setlocal path=.,**
    :autocmd FileType python setlocal wildignore=*/__pycache__/*,*.pyc
augroup end

"}}}

