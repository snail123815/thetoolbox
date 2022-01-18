set backspace=indent,eol,start
set history=50
set ruler
set showcmd
set incsearch
set hlsearch
set clipboard=unnamed
map Q gq
syntax on
filetype plugin indent on
set tabstop=4
set shiftwidth=4
set softtabstop=4
set expandtab
set autoindent
autocmd BufRead * autocmd FileType <buffer> ++once
      \ if &ft !~# 'commit\|rebase' && line("'\"") > 1 && line("'\"") <= line("$") | exe 'normal! g`"' | endif


command DiffOrig vert new | set bt=nofile | r ++edit # | 0d_ | diffthis
		  \ | wincmd p | diffthis

:set number

:augroup numbertoggle
:  autocmd!
:  autocmd BufEnter,FocusGained,InsertLeave,WinEnter * if &nu && mode() != "i" | set rnu   | endif
:  autocmd BufLeave,FocusLost,InsertEnter,WinLeave   * if &nu                  | set nornu | endif
:augroup END


set nolangremap
set backup
set backupdir-=.
set backupdir^=~/.nvimbk


let data_dir = has('nvim') ? stdpath('data') . '/site' : '~/.vim'
if empty(glob(data_dir . '/autoload/plug.vim'))
  silent execute '!curl -fLo '.data_dir.'/autoload/plug.vim --create-dirs  https://raw.githubusercontent.com/junegunn/vim-plug/master/plug.vim'
  autocmd VimEnter * PlugInstall --sync | source $MYVIMRC
endif


call plug#begin(has('nvim') ? stdpath('data') . '/plugged' : '~/.local/neovim/plugged')

" Declare the list of plugins.
"Plug 'davidhalter/jedi-vim'
Plug 'tpope/vim-sensible'
Plug 'ncm2/ncm2'
Plug 'ncm2/ncm2-jedi'
Plug 'roxma/nvim-yarp'
Plug 'ncm2/ncm2-bufword'
Plug 'ncm2/ncm2-path'
Plug 'preservim/tagbar'
Plug 'junegunn/seoul256.vim'

Plug 'mzlogin/vim-markdown-toc'

Plug 'godlygeek/tabular'
Plug 'plasticboy/vim-markdown'


" List ends here. Plugins become visible to Vim after this call.
call plug#end()

" enable ncm2 for all buffers
autocmd BufEnter * call ncm2#enable_for_buffer()

" IMPORTANT: :help Ncm2PopupOpen for more information
set completeopt=noinsert,menuone,noselect

" CTRL-C doesn't trigger the InsertLeave autocmd . map to <ESC> instead.
inoremap <c-c> <ESC>

" When the <Enter> key is pressed while the popup menu is visible, it only
" hides the menu. Use this mapping to close the menu and also start a new
" line.
inoremap <expr> <CR> (pumvisible() ? "\<c-y>\<cr>" : "\<CR>")

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

" Use <TAB> to select the popup menu:
inoremap <expr> <Tab> pumvisible() ? "\<C-n>" : "\<Tab>"
inoremap <expr> <S-Tab> pumvisible() ? "\<C-p>" : "\<S-Tab>"

let g:python3_host_prog = 'python3'

" 'plasticboy/vim-markdown'
let g:vim_markdown_toc_autofit = 1
let g:vim_markdown_math = 1

" 'mzlogin/vim-markdown-toc'
let g:vmt_dont_insert_fence = 1
let g:vmt_list_item_char = "-"


let g:netrw_banner = 0
let g:netrw_browse_split = 4
" 1 - open files in a new horizontal split
" 2 - open files in a new vertical split
" 3 - open files in a new tab
" 4 - open in previous window
let g:netrw_winsize = 80
let g:netrw_altv = 1
let g:netrw_liststyleW = 1

" clipboard
nnoremap Y "+y
vnoremap Y "+y
nnoremap yY ^"+y$

augroup filetype_md
	autocmd!
	:autocmd FileType markdown setlocal tabstop=2 shiftwidth=2 softtabstop=2 expandtab autoindent
	:autocmd FileType markdown noremap <silent> k gk
	:autocmd FileType markdown noremap <silent> j gj
	:autocmd FileType markdown noremap <silent> 0 g0
	:autocmd FileType markdown noremap <silent> $ g$
augroup end


inoremap <expr> <CR> InsertMapForEnter()
function! InsertMapForEnter()
    if pumvisible()
        return "\<C-Y>"
    elseif strcharpart(getline('.'),getpos('.')[2]-1,1) == '}'
        return "\<CR>\<Esc>O\<tab>123"
    elseif strcharpart(getline('.'),getpos('.')[2]-1,2) == '</'
        return "\<CR>\<Esc>O\<tab>abc"
    else
        return "\<CR>"
    endif
endfunction
