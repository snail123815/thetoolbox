#! /bin/zsh

# durand.dc@qq.com
# By Chao DU (杜超)
##############################################

# give pathes to different versions of inkscape
# change accrodingly if other wise
inkscape0_9="/Users/duc/Applications/Inkscape.app/Contents/Resources/bin/inkscape"
inkscape1_0="/usr/local/bin/inkscape"

# default values
version=1
format=jpg
dpi=300
quality=50

# help string, synopsis
showHelp() {
cat << EOF
Usage: ${1##*/} [options] FILE/DIR
Convert input FILE or the .svg files in DIR to .png or .jpg

inkscape and will be used to convert to .png
imagemagick is needed to convert the resulted .png to .jpg

Options:
    -v|--version
        inkscape version to use. Default $version
            0.9 = $inkscape0_9
            1   = $inkscape1_0
    -f|--fmt
        output format, png, jpg, pdf. Default $format
    --dpi
        resolution of png and jpg. Default 300
    -q|--quality
        quality setup of png to jpg. Default 50
EOF
}

convertSvgPng() {
    # $1 : absolute filepath
    # inkscape0.92 seems to take only absolute path
    # don't know about version 1.0, doesn't hurt
    if [[ ! $#=2 ]]; then
        echo "Exactly 2 arguments allowed for convertSvgPng()"
        exit 1
    fi
    local abssvg=$1
    abspng="${abssvg%.*}.png"
    if [[ ${version} = 0.9 ]]; then
        # for the propose of eval double quotes, the quotes needs to be escaped nicely
        local cmd="$inkscape0_9 -z -f \"$abssvg\" -e \"$abspng\" -d $dpi -C -b white"
    else
        local cmd="$inkscape1_0 -o \"$abspng\" -d $dpi -C -b white -y 1 \"$abssvg\""
    fi
    echo; echo "===================================================================="
    echo $cmd; echo
    eval "$cmd" # eval needs double quotes to work properly with such string
}

convertPngJpg() {
    # $1 : png filepath
    echo; echo "--------------------------------"; echo
    local cmd="magick \"$1\" -quality $quality \"${1%.*}.jpg\" && rm \"$1\""
    echo $cmd; echo
    eval "$cmd"
}

convertSvgPdf() {
    # $1 : absolute filepath
    # inkscape0.92 seems to take only absolute path
    # don't know about version 1.0, doesn't hurt
    if [[ ! $#=2 ]]; then
        echo "Exactly 2 arguments allowed for convertSvgPng()"
        exit 1
    fi
    local abssvg=$1
    abspdf="${abssvg%.*}.pdf"
    if [[ ${version} = 0.9 ]]; then
        # for the propose of eval double quotes, the quotes needs to be escaped nicely
        local cmd="$inkscape0_9 -z -f \"$abssvg\" -A \"$abspdf\" --export-text-to-path --export-ignore-filters -C -b white"
    else
        local cmd="$inkscape1_0 -o \"$abspdf\" -d $dpi -C -b white -y 1 --export-text-to-path --export-ignore-filters \"$abssvg\""
    fi
    echo; echo "===================================================================="
    echo $cmd; echo
    eval "$cmd" # eval needs double quotes to work properly with such string
}

getAbsFilePath() {
    # $1 : file path, doesn't matter relative or absolute
    echo $(cd "$(dirname "$1")" && pwd)/$(basename "$1")
}

# Start parsing options
while :; do
    case $1 in # starting parse the first switch
        -h|-\?|--help)
            showHelp $0 # Display a usage synopsis.
            exit 0
            ;;
        -v|--version) # Takes an option argument; ensure it has been specified.
            # now $1 is -v|--version
            version=$2
            # test string in a list, note the direction
            if [[ ! " 0.9 1 " = *" "$version" "* ]]; then
                echo "Wrong version number. 0.9 and 1 accepted" && exit 1
            fi
            shift # remove $1, now $1=VERSION, VERY IMPORTANT
            ;;
        -f|--fmt)
            # now $1 is -f|--fmt
            format=$2
            if [[ ! " png jpg pdf " = *" "$format" "* ]]; then # like the `in` operator
                echo "Wrong format, only png, jpg, pdf accepted" && exit 1
            fi
            shift # remove $1, now $1=FMT, VERY IMPORTANT
            ;;
        -q|--quality)
            # test int numbers only
            if (( $2 < 1 || $2 > 100 )); then
                echo "Wrong quality number $2, accepted range 1-100"; exit 1
            fi
            quality=$2
            shift
            ;;
        --dpi)
            # test if it is a number, float accepted:
            if ! [[ $2 =~ ^[0-9]+\.?[0-9]*$ ]]; then 
                echo "not a number $2"; exit 1
            # test if in or out of a range, use bc:
            elif [[ $(bc <<<"($2 < 0.1) || ($2 > 10000)") = 1 ]]; then 
                echo "dpi range 0.1-10000"; exit 1
            fi
            dpi=$2
            shift
            ;;
        --) # specified end of all options if given.
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            # put the warning into stderr (&2)
            shift
            ;;
        *)  
            # If it sees the same thing as previos loop,
            # that means no more options, so break out of the loop.
            break
    esac
    # IMPORTANT:
    shift # if case is iterating to an argument value, shift and case again
done

# parse actual file or dir from left over argument
pathInput=$(getAbsFilePath "$@")

# generate an array of files to process
files=()
for var in "$@"; do
    pathInput=$(getAbsFilePath "$var")
    if [ -d $pathInput ]; then
        for fsvg in "${pathInput}/"*.svg; do
            files+=$fsvg
        done
    elif [ -f $pathInput ]; then
        files+=$pathInput
    fi
done

# process
for fsvg in ${files[@]}; do
    if [[ $format = pdf ]]; then
        convertSvgPdf "$fsvg"
    else # png or jpg
        convertSvgPng "$fsvg"
        if [[ $format = jpg ]]; then
            convertPngJpg "$abspng"
        fi
    fi
done