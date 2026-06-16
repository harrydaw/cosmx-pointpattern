param(
    [Parameter(Mandatory=$true)][string]$Pptx,
    [Parameter(Mandatory=$true)][string]$OutPng,
    [int]$Width = 2480
)
$Pptx = (Resolve-Path $Pptx).Path
$height = [int]($Width * 118.9 / 84.1)
$ppt = New-Object -ComObject PowerPoint.Application
try {
    $pres = $ppt.Presentations.Open($Pptx, $true, $false, $false)  # ReadOnly, Untitled, WithWindow=false
    $pres.Slides.Item(1).Export($OutPng, "PNG", $Width, $height)
    $pres.Close()
    "Exported $OutPng ($Width x $height)"
} finally {
    $ppt.Quit()
    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($ppt) | Out-Null
}
