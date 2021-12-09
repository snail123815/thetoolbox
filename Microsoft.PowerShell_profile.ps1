New-Alias ll Get-ChildItem

Remove-Alias r # for removing the system R and use only conda R if that is the purpose

function rsyncps([string]$from, [string]$to, [string]$wsl='ubuntu-wsl1')
{
    # Purpose is to have a restartable copy function. robocopy can do it but the
    # transfer behaviour is very weird. Without specifying -MT > 32, it won't use full
    # bandwidth. Also, if you specify -MT in any number, the disk will be fully occupied
    # by it. Which makes no sense.

    # pass in a wsl version 1 to ensure efficiency (probably don't matter)
    # you need to mount all accessible drives in /mnt/ before using this function, or 
    # we won't find target/source location
    # First make dir for mount point, then add lines in /etc/fstab using sudo
    # J: /mnt/j/ drvfs defaults 0 0
    # P: /mnt/p/ drvfs defaults 0 0
    # this should be permenent if your windows version >= 17093
    if ($from -match '^[a-z]\:[\\\/]') {
        $from = $from -replace '^[a-z]\:[\\\/]', "/mnt/$($Matches[0][0])/".ToLower()
    }
    if ($to -match '^[a-z]\:[\\\/]') {
        $to = $to -replace '^[a-z]\:[\\\/]', "/mnt/$($Matches[0][0])/".ToLower()
    }
    $from = $from -replace '\\', "/"
    $from = $from -replace '^\.\/',""
    $to = $to -replace '\\', "/"
    $to = $to -replace '^\.\/',""
    $cmd = "wsl -d $wsl rsync -a --info=progress2 `"$from`" `"$to`""
    Write-Output $cmd
    Invoke-Expression $cmd
}

function vscode ($path="")
{
    $code_path="`"${Env:LOCALAPPDATA}\Programs\Microsoft VS Code\Code.exe`""
    $cmd = "Start-Process -FilePath $code_path"
    if ($path -ne "") {
        $cmd = $cmd + " -ArgumentList `"```"$path```"`""
    }
    Write-Output $cmd
    Invoke-Expression $cmd
}

function tarpigz ($path)
{
    $path = $path -replace '^\.\\',""
    # If it is the local path, powershell often add .\ in the beginning, this will
    # not be recognised by tar in linux probably because the whole thing is not
    # passed to bash but just as an argument.
    $target = $path + '.tar.gz'
    wsl -d ubuntu-wsl1 tar cvf $target --use-compress-program=pigz $path
    # do not pipe because that will pipe back to powershell
    # maybe try to use bash and "tar ..." and see if that can parse the pip correctly
}