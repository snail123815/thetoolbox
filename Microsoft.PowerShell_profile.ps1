New-Alias ll Get-ChildItem

Remove-Alias r # for removing the system R and use only conda R if that is the purpose

function rsyncps() {
    param(
        [Parameter(Mandatory=$False)][string]$wsl='ubuntu-wsl1', # wsl distribution name
        [Parameter(
            Mandatory=$True,
            Position=1, # To make sure wsl is only set when specifying -wsl
            ValueFromRemainingArguments=$True
        )][string[]]$paths
    )
    # Purpose is to have a restartable copy function. robocopy can do it but the
    # transfer behaviour is very weird. Without specifying -MT > 32, it won't use full
    # bandwidth. Also, if you specify -MT in any number, the disk will be fully occupied
    # by it. Which makes no sense.
    # rysnc is better in utilising all bandwidth at all times. However it is not perfect.
    # For some reason it uses bandwidth alternatively with
    # 'Windows Defender Advanced Threat Protection Service', occupancy is half-half.
    # Albeit that, disk read and network sending looks fine most of the time, and speed is
    # up-to-expectation: Showing speed ~12 MB/s + gaps in transferring buffer ~= 8.4 MB/s on a
    # 100 Mbps (12.5 MB/s max) band. (good enough for both transfer and check)

    # pass in a wsl version 1 to ensure efficiency (probably don't matter)
    # you need to mount all accessible drives in /mnt/ before using this function, or
    # we won't find target/source location
    # First make dirs for mount point, then add lines in /etc/fstab using sudo
    # J: /mnt/j/ drvfs defaults 0 0
    # P: /mnt/p/ drvfs defaults 0 0
    # this should be permenent if your windows version >= 17093
    $cmd = "wsl -d $wsl "
    $froms = $paths[0..($paths.Length - 2)]
    $to = $paths[-1]
    for ($i = 0; $i -le ($froms.length-1); $i += 1) {
        if ($froms[$i] -match '^[a-z]\:[\\\/]') {
            $froms[$i] = $froms[$i] -replace '^[a-z]\:[\\\/]', "/mnt/$($Matches[0][0])/".ToLower() `

        }
        $froms[$i] = $froms[$i] -replace '\\', "/" -replace '^\.\/',""
    }
    if ($to -match '^[a-z]\:[\\\/]') {
        #$netdrives = 'J:', 'K:'
        #$netdrives = ((net use | select-object -skip 3) -replace '\s{2,}', ' ' -replace '-', ''
        #             | ConvertFrom-Csv -delimiter ' '
        #             | Where-Object {$_.Local -like '*:'}).Local
        # sudo is needed if target location is a network location (mounted to windows)
        # another solution to get netdrives neatly but lags a bit:
        $netdrives = (Get-CimInstance -Class Win32_LogicalDisk | Where-Object {$_.DriveType -eq 4}).DeviceID
        if ($to.Substring(0,2) -in $netdrives) { $cmd += 'sudo ' }
        $to = $to -replace '^[a-z]\:[\\\/]', "/mnt/$($Matches[0][0])/".ToLower()
    }
    $to = $to -replace '\\', "/" -replace '^\.\/',""
    $cmd += "rsync -a --info=progress2"
    foreach ($from in $froms) {
        $cmd += " `"$from`""
    }
    $cmd += " `"$to`""
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
    $cmd = "wsl -d ubuntu-wsl1 tar cvf `"$target`" --use-compress-program=pigz `"$path`""
    Write-Output $cmd
    Invoke-Expression $cmd
    # maybe try to use bash and "tar ..." and see if that can parse the pip correctly
}