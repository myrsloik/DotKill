# DotKill
Spatial dotcrawl remover for VapourSynth

# Functions
dotkill.DotKill(clip clip, int napply=1, bint ignorematch=0)

napply: The number of times to apply the internal filter. Usally a number between 1 and 4 will have the best results and using too high values may cause artifacting.

ignorematch: If true then matching hints from VFM are ignored when processing. This argument mostly exists for debugging purposes and shouldn't be changed.

# Usage
DotKill can be used after fieldmatching and decimation unlike many other filters. Note that currently only VFM provides the hints so placing it after other field matchers will fail.

clip = core.vivtc.VFM(clip)

clip = core.vivtc.VDecimate(clip)

clip = core.dotkill.DotKill(clip, napply=4)

If you don't need IVTC then simply use DotKill by itself.
