notebookSync() {
	echo
	local CWD=$PWD
	local NOTEBOOK=~/jianguoYUN/Reports/NOTEBOOK/
	echo "PATH to notebook:"
	echo $NOTEBOOK
	echo
	local COMMENT=${1:-$(date)}
	echo Comment is: "$COMMENT"
	if [ "$COMMENT" = "" ]; then
		echo succeed
	fi
	cd $NOTEBOOK
	git pull
	git add .
	git commit -m $COMMENT
	git push
	cd $CWD
	echo
    echo "####### Sync images ####"
    echo
    imgbedSync
    echo
	echo "####### ALL DONE #######"
	echo
}

imgbedSync() {
	local currentdir=$(PWD)
    local imgbedRepo="$HOME/GitProjects/imgRepo"
    echo "Images local path: "$imgbedRepo
    echo
	cd $imgbedRepo && git pull -v
    echo
    #ssh-add ~/.ssh/gitee_id_rsa_misc
    git push -v ge
	cd $currentdir
}

jdk() {
        local version=$1
        export JAVA_HOME=$(/usr/libexec/java_home -v "$version");
        echo "Current JDK version:"
        java -version
	echo
	echo "Avaliable JDKs:"
	/usr/libexec/java_home -V
 }

convertPics() {
    local resize=0
    local reformat=0
    local sourceExt=0
    local RM=0
    local args=()
    local moreArgs=()
    local showHelp() {
        cat << EOF
Usage: convertPics [options] FILE(s)/DIR(s)
Convert input img file or files in directory into target format.

inkscape will be used to convert .svg files.

Options:
    Key value:
    -rf     TARGET_FORMAT_EXTENSION
            If you want to convert format. Will do "magick convert".
            Eg. "jpg"; "png"
    -rs     RESIZE_FACTOR
            If you want to resize picutre. Will add "-resize RF" in
            "magick convert" command. Eg. "40%"; "1000"(width); "x1000"(height);
    -sf     SOURCE_FORMAT_EXTENSION
            When giving DIR as input, this parameter specifies images of which
            format will be converted. Mandatory when giving DIR input
    
    When converting .svg files:
    -dpi    DPIvalue
            Default to 300
    -t|-trans|-transparent
            A switch, default is off. When converting to PNG, whether to keep transparency.
    
    Other switches:
    -c      Remove source file
    -cy     Remove source file without asking
    -nr     Non-recursive, only valid when passing DIR with sub directories. Default is recursive.
    -h      Show this help message and exit

    Other unreconised parameters will be passed to inkscape or magick

EOF
    }

    while [ ! -f $1 ] && [ ! -d $1 ]; do
        case $1 in
            -h)
                showHelp
                return ;;
            -rf)
                reformat=1
                local targetExt=$2
                args+=("-rf" "$targetExt")
                shift 2 ;;
            -rs)
                resize=1
                local resizeFactor=$2
                args+=("-rs" "$resizeFactor")
                shift 2 ;;
            -sf)
                sourceExt=$2
                args+=("-sf" "$sourceExt")
                shift 2 ;;
            -c)
                echo "With -c, the original file(s) will be removed, y to confirm."
                echo "Press other keys to ignore, files will not be removed."
                if [[ $0 == *'zsh' ]]; then
                    read -q CONTINUE # for zsh
                else
                    read CONTINUE # for BASH, I don't know others...
                fi
                echo # new line
                if [ "$CONTINUE" = "y" ]; then
                    RM=1
                    args+="-cy"
                fi
                shift 1 ;;
            -cy)
                RM=1
                args+="-cy"
                shift 1 ;;
            *)
                moreArgs+=$1
                shift 1 ;;
        esac
    done

	for file in $@; do
        if [ -f $file ]; then
            local name="${file%.*}"
            local ext="${file##*.}"
            if [ ! $sourceExt = "$ext" ] && [ ! $sourceExt = 0 ]; then
                continue
            fi

            if [ $ext = 'svg' ]; then
                local CMD="convertSvg "
                if [ ! $reformat = 1 ]; then
                    echo "input svg without reformat will anyway be converted to png"
                    reformat=1
                    local targetExt="png"
                    args+=("-rf" "$targetExt")
                fi

                if [[ "${#args[@]}" > 0 ]]; then
                    CMD+=" ${args[@]}"
                fi
                if [[ "${#moreArgs[@]}" > 0 ]]; then
                    CMD+=" ${moreArgs[@]}"
                fi
                CMD+=" \"$file\""
            else
                local CMD="magick convert"
                if [[ "${#moreArgs[@]}" > 0 ]]; then
                    CMD+=" ${moreArgs[@]}"
                fi
                CMD+=" \"$file\""

                local newFile=${file%.*}.$ext
                if [ $resize = 1 ]; then
                    CMD+=" -resize"" $resizeFactor"
                    newFile=${newFile%.*}.resize.$ext
                fi
                if [ $reformat = 1 ]; then
                    if [ $targetExt = $ext ]; then
                        echo "reformat option error, target extension same as input"
                        exit 1
                    fi
                    newFile=${newFile%.*}.$targetExt
                fi
                CMD+=" \"$newFile\""
            fi

            if [[ $CMD == "magick"* ]]; then
                echo; echo $CMD; echo
            fi
            eval "$CMD"

            if [ $RM = 1 ]; then
                rm $file
            fi
        elif [ -d $file ]; then
            if [ $sourceExt = 0 ]; then
                echo "Please provide \"-sf SOURCE_FORMAT_EXTENSION\" as arguments when giving DIR"
                return
            fi
            local CMD="_convertPicsInDirs ${args[@]} ${moreArgs[@]} \"$file\""
            #echo $CMD
            eval "$CMD"
        fi

    done
}


convertSvg() {
    local inkscapePath="/usr/local/bin/inkscape"
    local dpi=300
    local targetExt="png"
    local whiteBg=1
    local moreArgs=()
    
    while [ ! -f $1 ]; do
        case $1 in
            -dpi)
                dpi=$2
                shift 2 ;;
            -rf) # reformat
                targetExt=$2
                shift 2 ;;
            -transparent | -trans | -t)
                whiteBg=0
                shift 1 ;;
            -sf) # ignore
                shift 2 ;;
            *)
                moreArgs+=$1
                shift 1 ;;
        esac
    done

    local svgPath=$(_getAbsFilePath "$1")
    
    local targetPath="${svgPath%.*}.$targetExt"

    local CMD="$inkscapePath -o \"$targetPath\" -d $dpi -C"
    if [ $whiteBg = 1 ]; then
        CMD+=" -b white"
    fi
    CMD+=" \"$svgPath\""
    echo; echo $CMD; echo
    eval "$CMD" # eval needs double quotes to work properly with such string
    if [ ! $targetExt = "png" ] && [ ! $targetExt = "pdf" ]; then
        CMD="convertPics -rf $targetExt -cy ${moreArgs[@]} \"$targetPath\""
        #echo $CMD
        eval "$CMD"
    fi
}

compressPDF() {
	local sourceFile=$1
	local targetFile=${1%.*}.resize.pdf
	gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer -sOutputFile=$targetFile $sourceFile
}


_getAbsFilePath() {
    # $1 : file path, doesn't matter relative or absolute
    echo $(cd "$(dirname "$1")" && pwd)/$(basename "$1")
}


_convertInOneDir() {
    local moreArgs=()
    while [ ! -d $1 ]; do
        moreArgs+=$1
        shift 1
    done
    for file in $1/*.*; do
        local CMD="convertPics ${moreArgs[@]} \"$file\""
        #echo $CMD
        eval "$CMD"
    done
}


_convertPicsInDirs() {
	local moreArgs=()
	local RECURSIVE=1
	if [ $1 = "-nr" ]; then
		RECURSIVE=0
		shift 1
	fi
	local RM=0
	while [ ! -d $1 ]; do
		moreArgs+=$1
		shift 1
	done

	if [ $RECURSIVE = 0 ]; then
        local CMD="_convertInOneDir ${moreArgs[@]} \"$1\""
        #echo $CMD
        eval "$CMD"
	else
		local OIFS=$IFS
		local IFS=$'\n' DIRS=($(find "$1" -type d))
		local IFS=$OIFS
		for d in ${DIRS}; do
            local CMD="_convertInOneDir ${moreArgs[@]} \"$d\""
            #echo $CMD
            eval "$CMD"
		done
	fi
}


