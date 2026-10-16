param(
    [string]$Binary = "$PSScriptRoot/../build/benchmarks/Release/bench_2d_compare.exe",
    [string]$OutputDirectory = "$PSScriptRoot/../build/comparison",
    [int]$Frames = 300,
    [int]$Count = 1000
)
$ErrorActionPreference = 'Stop'
$null = New-Item -ItemType Directory -Force -Path $OutputDirectory
# Run serially, after builds/tests finish. Each invocation alternates engine
# order over three rounds. Baseline means current corrected code with selected
# optimization flags disabled, not an older Butter revision.
$runs = @(
    @('all', 8, 1, 'awake', 'baseline'),
    @('all', 8, 1, 'awake', 'optimized'),
    @('all', 8, 1, 'natural', 'optimized'),
    @('stack', 4, 4, 'awake', 'optimized'),
    @('mass_ratio', 4, 4, 'awake', 'optimized'),
    @('stack', 4, 4, 'awake', 'soft'),
    @('mass_ratio', 4, 4, 'awake', 'soft'),
    @('stack', 8, 1, 'awake', 'anchors'),
    @('circles', 8, 1, 'awake', 'simd'),
    @('circles', 8, 1, 'awake', 'parallel')
)
foreach ($run in $runs) {
    $name = $run -join '-'
    & $Binary $run[0] $Frames $Count $run[1] $run[2] $run[3] $run[4] `
        1> "$OutputDirectory/$name.csv" 2> "$OutputDirectory/$name.metadata.txt"
    if ($LASTEXITCODE -ne 0) { throw "Comparison failed: $name" }
    Write-Output "Completed $name"
}
