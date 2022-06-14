# af-magic.zsh-theme
# Repo: https://github.com/andyfleming/oh-my-zsh
# Direct Link: https://github.com/andyfleming/oh-my-zsh/blob/master/themes/af-magic.zsh-theme


# settings
typeset +H return_code="%(?..%{$fg[red]%}╰%? %{$reset_color%})"
typeset +H my_gray="$FG[242]"
typeset +H my_orange="$FG[214]"
# separator dashes size
function afmagic_dashes {
	local PYTHON_ENV="$VIRTUAL_ENV"
    local HOST_NAME=$(hostname | cut -f 1 -d .)
    local CURRENT_PATH="${PWD/#$HOME/'~'}"
    local CURRENT_GIT="$(git_prompt_info)$(hg_prompt_info) "
    local zero='%([BSUbfksu]|([FK]|){*})' # https://stackoverflow.com/a/10564427/6823079
	[[ -z "$PYTHON_ENV" ]] && PYTHON_ENV="$CONDA_DEFAULT_ENV"

	if [[ -n "$PYTHON_ENV" && "$PS1" = *\(* ]]; then
		echo $(( COLUMNS - ${#PYTHON_ENV} - 5 - ${#USER} - ${#HOST_NAME} - ${#CURRENT_PATH} - ${#${(S%%)CURRENT_GIT//$~zero/}} ))
	else
		echo $(( COLUMNS - 2 - ${#USER} - ${#HOST_NAME} - ${#CURRENT_PATH} - ${#${(S%%)CURRENT_GIT//$~zero/}} ))
	fi
}

# primary prompt
#PS1='$my_gray%n@%m%{$reset_color%} $FG[032]%~$(git_prompt_info)$(hg_prompt_info) $FG[237]${(l.$(afmagic_dashes)..-.)}%{$reset_color%}
PS1='$my_gray%n@%m%{$reset_color%} $FG[032]%~$(git_prompt_info)$(hg_prompt_info) $FG[237]${return_code%}%{$reset_color%}
$FG[105]%(!.#.$)%{$reset_color%} '
#RPS1='${return_code}'

# git settings
ZSH_THEME_GIT_PROMPT_PREFIX="$FG[075]($FG[078]"
ZSH_THEME_GIT_PROMPT_CLEAN=""
ZSH_THEME_GIT_PROMPT_DIRTY="$my_orange*%{$reset_color%}"
ZSH_THEME_GIT_PROMPT_SUFFIX="$FG[075])%{$reset_color%}"

# hg settings
ZSH_THEME_HG_PROMPT_PREFIX="$FG[075]($FG[078]"
ZSH_THEME_HG_PROMPT_CLEAN=""
ZSH_THEME_HG_PROMPT_DIRTY="$my_orange*%{$reset_color%}"
ZSH_THEME_HG_PROMPT_SUFFIX="$FG[075])%{$reset_color%}"

# virtualenv settings
ZSH_THEME_VIRTUALENV_PREFIX=" $FG[075]["
ZSH_THEME_VIRTUALENV_SUFFIX="]%{$reset_color%}"
