
notebookSync() {
	local CWD=$PWD
    local NOTEBOOKs=(
        ~/NOTEBOOK/
        ~/GitProjects/IBLProteomics
    )
    for n in ${NOTEBOOKs[@]}; do
        _notebookSync $n
    done
	cd $CWD

	echo
    echo "####### Sync images ####"
    echo
    imgbedSync
    echo
	echo "####### ALL DONE #######"
	echo
}

_notebookSync() {
	echo
	local NOTEBOOK=$1
	echo "PATH to notebook:"
	echo $NOTEBOOK
	echo
	local COMMENT=${:-$(date)}
	echo Comment is: "$COMMENT"
	if [ "$COMMENT" = "" ]; then
		echo succeed
	fi
	cd $NOTEBOOK
	git pull
	git add .
	git commit -m $COMMENT
	git push
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

cpics() {
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
                local CMD="csvg "
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


csvg() {
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
    
    if [ ! $targetExt = "png" ] && [ ! $targetExt = "pdf" ]; then
        local targetPath="${svgPath%.*}.png"
    else
        local targetPath="${svgPath%.*}.$targetExt"
    fi

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


cpdf() {
    "
    # Usage:
    # cpdf input.pdf
    #     will generate input.resize.pdf with dpi=300 quality=printer
    # cpdf -png input.pdf
    #     will generate input-N.png file for each page and then combine them into input.png.pdf file
    # cpdf -jpeg/jpg input.pdf
    #     will generate input-N.jpeg file for each page and then combine them into input.jpeg.pdf file
    #
    # you can add -r 300 ro -res 300 when using -png or -jpeg/jpg to specify resolution of generated picture
    # Requirements: (1) Ghostscript needs to be installed on the local system.
    #              (2) ImageMagick needs to be installed on the local system.
    #
    "
    local sDEVICE=pdfwrite
    local ext=pdf
    local CMD="gs -dNOPAUSE -dBATCH"
    local res=" -r600"
    while [ ! -f $1 ]; do
        case $1 in
            -png)
                sDEVICE=png16m
                ext=png
                shift 1 ;;
            -jpg | -jpeg)
                sDEVICE=jpeg
                ext=${1:1}
                shift 1 ;;
            -r | -res)
                res=" -r"$2
                shift 2 ;;
            *)
                if [ ! -f $1 ]; then
                    echo "File "\"$1\"" not found."
                    return 1
                fi
        esac
    done

    local sourceFile=$1
    local sourceExt="${sourceFile##*.}"
    if [ ! $sourceExt = 'pdf' ]; then
        echo "Target file should be *.pdf"
        return 1
    fi

    CMD+=" -sDEVICE=$sDEVICE"

    if [ $ext = 'pdf' ]; then
        CMD+=" -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer"
        local targetFile=${sourceFile%.*}.resize.$ext
    else
        CMD+=$res
        local targetFile=${sourceFile%.*}-%03d.$ext
    fi

    CMD+=" -sOutputFile=""\"$targetFile\""" ""\"$sourceFile\""

    echo $CMD
    eval "$CMD"

    if [ ! $ext = 'pdf' ]; then
        local targetPdf=${sourceFile%.*}.$ext.pdf
        local $I2PCMD="magick convert"
        I2PCMD+=" \"${sourceFile%.*}*.$ext\""
        I2PCMD+=" ""\"$targetPdf\""
        echo $I2PCMD
        eval "$I2PCMD"
    fi

}


__join_by() {
    local d=$1 f=$2
    printf $s $f
    if shift 2; then
        for i in $@; do
            printf %s "$d""$i"
        done
    fi
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


