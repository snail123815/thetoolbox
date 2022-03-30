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

resizePics() {
    # resize pictures as arguments
    # Usage:
    # resizePics 20% [-c] [resize args] img.jpg img1.jpg img2.jpg ...
    # -c - clean up after resize. DANGEROUS!
	local resize=$1
	shift 1

	local REMOVEORIGINAL=0

	if [ $1 = "-c" ]; then
		REMOVEORIGINAL=1
		shift 1
	fi
	local moreArgs=()
	while [ ! -f $1 ]; do
		moreArgs+=$1
		shift 1
	done

	for file in $@; do
		local name="${file%.*}"
		local ext="${file##*.}"
		local newFile="${name}".resize."${ext}"
		echo magick convert $file -resize $resize ${moreArgs[@]} $newFile
		magick convert $file -resize $resize ${moreArgs[@]} $newFile
		if [ $REMOVEORIGINAL = 1 ]; then
			echo "With -c, the original file $file will be removed, y to continue."
			#read CONTINUE # for BASH
			read -q CONTINUE # for zsh
			echo # for zsh
			if [ "$CONTINUE" != "y" ]; then
				return
			else
				rm $file
			fi
		fi
	done
}

resizePicInDir() {
    # Resize images in a folder
    # Usage:
    # resizePicInDir jpg png [-nr] [-c] [resize args] ./target/dir/
    # source file extension as 1st argument
	# target file extension as 2nd argument
    # -nr recursive to subfolders, defaults to NOT recur
    # -c clean up after resize. DANGEROUS! will remove source file.

    # The actual function for one dir
	resizeOnedir() {
        # Usage:
        # resizeOnedir $REMOVEORIGINAL ./path/
        # #REMOVEORIGINAL is either 1 or 0

		local RM=$1 # Receive $REMOVEORIGINAL
		shift 1

        local moreArgs=()
        while [ ! -d $1 ]; do
            moreArgs+=$1
            shift 1
        done

		# See if file with specified extension exists
		X=$(find $1 -maxdepth 1 -name "*.$format")

		# $X is just a string, not an array.
		# Converting it into an array is painful.
		# Read -a only works with BASH
		# ${=X} only work in ZSH, but need to change IFS beforehand.
        # IFS is a dangerous argument. Try to avoid it...
		
		if [ ! -z $X ]; then
			for file in $1/*.$format; do
				if [ -d $file ]; then
					continue
				fi
				newFile=${file%.*}.resize."$format"
                CMD = "magick convert $file -resize $resize ${moreArgs[@]} $newFile"
                eval "$CMD"
				if [ $RM = 1 ]; then
					echo rm $file
					rm $file
				fi
			done
		fi
	}

    # Now starts the main function
	local resize=$1
	shift 1
	local format=$1
	shift 1
	local moreArgs=()
	local RECURSIVE=0
	if [ ! $1 = "-nr" ]; then
		RECURSIVE=1
		shift 1
	fi
	local REMOVEORIGINAL=0
	if [ $1 = "-c" ]; then
		REMOVEORIGINAL=1
		echo "With -c, the original file will be removed, y to continue."
        if [[ $0 == *'zsh' ]]; then
            read -q CONTINUE # for zsh
        else
            read CONTINUE # for BASH
        fi
		echo # new line
		if [ "$CONTINUE" != "y" ]; then
			return
		fi
		shift 1
	fi
	while [ ! -d $1 ]; do
		moreArgs+=$1
		shift 1
	done

	if [ $RECURSIVE = 0 ]; then
		resizeOnedir $REMOVEORIGINAL ${moreArgs[@]} $1
	else
		local OIFS=$IFS
		local IFS=$'\n' DIRS=($(find "$1" -type d))
		IFS=$OIFS
		# Do forget about ZSH/BASH compatibility aaaaaah, also issues with IFS aaaaaah
		for d in ${DIRS}; do
			resizeOnedir $REMOVEORIGINAL ${moreArgs[@]} $d
		done
	fi

}

# TODO: write a unified function for magick convert and resize
# ON one file, ON more files, ON one dir, ON multiple dirs

_convertPics() {
    # including resize and change format
    local resize=0
    local reformat=0
    local RM=0
    local moreArgs=()
    if [ $1 = 'rs' ]; then # resize
        shift 1
        resize=1
        local resizeFactor=$1
        shift 1
    fi
    if [ $1 = 'rf' ]; then # reformat
        shift 1
        reformat=1
        local targetExt=$1
        shift 1
    fi
    if [ $1 = "-c" ]; then
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
		fi
        shift 1
    fi
	while [ ! -f $1 ]; do
		moreArgs+=$1
		shift 1
	done
	for file in $@; do
		local name="${file%.*}"
		local ext="${file##*.}"

        if [ $ext = 'svg' ]; then
            local CMD=""
        else
            local CMD="magick convert"
            if [[ "${#moreArgs[@]}" > 0 ]]; then
                CMD+=" ${moreArgs[@]}"
            fi
            CMD+=" $file"
            if [ $resize = 1 ]; then
                CMD+=" -resize"" $resizeFactor"
            fi
            if [ $reformat = 1 ]; then
                local newFile=${file%.*}.format.$targetExt
            else
                local newFile=${file%.*}.format.$ext
            fi
            CMD+=" $newFile"
        fi

        echo $CMD
        eval "$CMD"

        if [ $RM = 1 ]; then
            rm $file
        fi
    done
}

_getAbsFilePath() {
    # $1 : file path, doesn't matter relative or absolute
    echo $(cd "$(dirname "$1")" && pwd)/$(basename "$1")
}
f
_convertSvgPng() {
    local inkscapePath="/usr/local/bin/inkscape"

    # default values
    local dpi=300
    local quality=90
    # $1 : absolute filepath
    if [[ ! $#=2 ]]; then
        echo "Exactly 2 arguments allowed for convertSvgPng()"
        exit 1
    fi
    local svgPath=$1
    pngPath="${svgPath%.*}.png"
    local CMD="$inkscapePath -o \"$pngPath\" -d $dpi -C -b white -y 1 \"$svgPath\""
    echo $CMD; echo
    eval "$CMD" # eval needs double quotes to work properly with such string
}

convertPngJpg() {
    # $1 : png filepath
    echo; echo "--------------------------------"; echo
    local CMD="magick \"$1\" -quality $quality \"${1%.*}.jpg\" && rm \"$1\""
    echo $CMD; echo
    eval "$CMD"
}

convertSvgPdf() {
    # $1 : absolute filepath
    # inkscape0.92 seems to take only absolute path
    # don't know about version 1.0, doesn't hurt
    if [[ ! $#=2 ]]; then
        echo "Exactly 2 arguments allowed for convertSvgPng()"
        exit 1
    fi
    local svgPath=$1
    abspdf="${svgPath%.*}.pdf"
    if [[ ${version} = 0.9 ]]; then
        # for the propose of eval double quotes, the quotes needs to be escaped nicely
        local CMD="$inkscape0_9 -z -f \"$svgPath\" -A \"$abspdf\" --export-text-to-path --export-ignore-filters -C -b white"
    else
        local CMD="$inkscapePath -o \"$abspdf\" -d $dpi -C -b white -y 1 --export-text-to-path --export-ignore-filters \"$svgPath\""
    fi
    echo; echo "===================================================================="
    echo $CMD; echo
    eval "$CMD" # eval needs double quotes to work properly with such string
}

getAbsFilePath() {
    # $1 : file path, doesn't matter relative or absolute
    echo $(cd "$(dirname "$1")" && pwd)/$(basename "$1")
}

convertPics() {
	local targetFormat=$1
	shift 1
	local REMOVEORIGINAL=0
	if [ $1 = "-c" ]; then
		REMOVEORIGINAL=1
		shift 1
	fi


	local moreArgs=()
	while [ ! -f $1 ]; do
		moreArgs+=$1
		shift 1
	done

	for file in $@; do
		local name="${file%.*}"
		local ext="${file##*.}"
		local newFile=${file%.*}.format.$targetFormat
		CMD="magick convert ${moreArgs[@]} $file $newFile"
        echo $CMD
        eval "$CMD"
		if [ $REMOVEORIGINAL = 1 ]; then
			echo "With -c, the original file $file will be removed, y to continue."
			#read CONTINUE # for BASH
			read -q CONTINUE # for zsh
			echo # for zsh
			if [ "$CONTINUE" != "y" ]; then
				return
			else
				rm $file
			fi
		fi
	done
}

convertPicsInDir() {
	local sourceFormat=$1
	shift 1
	local targetFormat=$1
	shift 1
	local moreArgs=()
	local RECURSIVE=1
	if [ $1 = "-nr" ]; then
		RECURSIVE=0
		shift 1
	fi
	local REMOVEORIGINAL=0
	if [ $1 = "-c" ]; then
		REMOVEORIGINAL=1
		echo "With -c, the original file will be removed, y to continue."
		#read CONTINUE # for BASH
		read -q CONTINUE # for zsh
		echo # for zsh
		if [ "$CONTINUE" != "y" ]; then
			return
		fi
		shift 1
	fi
	while [ ! -d $1 ]; do
		moreArgs+=$1
		shift 1
	done
	convertInOneDir() {
		local RM=$1
		shift 1
		local X=$(find $1 -maxdepth 1 -name "*.$sourceFormat")

		if [ ! -z $X ]; then
			for file in $1/*.$sourceFormat; do
				if [ -d $file ]; then
					continue
				fi
				local newFile=${file%.*}.format.$targetFormat
				echo magick convert ${moreArgs[@]} $file $newFile
				magick convert ${moreArgs[@]} $file $newFile
				if [ $RM = 1 ]; then
					echo rm $file
					rm $file
				fi
			done
		fi
	}


	if [ $RECURSIVE = 0 ]; then
		convertInOneDir $REMOVEORIGINAL $1
	else
		local OIFS=$IFS
		local IFS=$'\n' DIRS=($(find "$1" -type d))
		local IFS=$OIFS
		for d in ${DIRS}; do
			convertInOneDir $REMOVEORIGINAL $d
		done
	fi
}

compressPDF() {
	local sourceFile=$1
	local targetFile=${1%.*}.resize.pdf
	gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer -sOutputFile=$targetFile $sourceFile
}
