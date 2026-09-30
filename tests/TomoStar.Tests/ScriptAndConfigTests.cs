// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Cli;
using TomoStar.Core.Model;

namespace TomoStar.Tests;

/// <summary>The configuration layers and the script interpreter.</summary>
public class ScriptAndConfigTests : IDisposable
{
    private readonly string _dir = Directory.CreateTempSubdirectory("tomostar-script-").FullName;

    public void Dispose()
    {
        try { Directory.Delete(_dir, true); } catch (IOException) { }
    }

    [Fact]
    public void ConfigurationLayersApplyInOrder()
    {
        var file = Path.Combine(_dir, "c.json");
        File.WriteAllText(file, """
            {
              // comments are allowed
              "Tomography": { "Smoothing": 7, "Iterations": 3 },
              "OutsideData": "WhenRaysCross",
            }
            """);
        var c = TomoConfig.Build(file, ["Tomography.Smoothing=11", "Attenuation.Phase=S", "LCurve.Dampings=[1,10]"], null);
        Assert.Equal(11, c.Tomography.Smoothing);       // --set wins over the file
        Assert.Equal(3, c.Tomography.Iterations);       // the file wins over the defaults
        Assert.Equal(20, c.Tomography.DampingVelocity); // the default
        Assert.Equal(OutsideData.WhenRaysCross, c.OutsideData);
        Assert.Equal(Phase.S, c.Attenuation.Phase);
        Assert.Equal([1.0, 10.0], c.LCurve.Dampings);
    }

    [Fact]
    public void TokensHonourQuotes()
    {
        Assert.Equal(["invert", "--name", "a b", "--set", "X=\"q\""], ScriptRunner.Tokenize("invert --name 'a b' --set \"X=\\\"q\\\"\""));
    }

    [Fact]
    public void ScriptsExpandVariablesLoopsConditionsAndJsonValues()
    {
        File.WriteAllText(Path.Combine(_dir, "values.json"), """{ "Recommended": { "Damping": 42.5 } }""");
        File.WriteAllText(Path.Combine(_dir, "run.tomo"), """
            # a script
            default OUT = out
            set D = ${json:values.json:Recommended.Damping}
            foreach M in ak135 iasp91
              model1d ${M} --out ${OUT}/${M}_d${D}.txt
            end
            if missing ${OUT}/never
              model1d PREM \
                --out ${OUT}/prem.txt
            end
            if exists ${OUT}/never
              model1d nonexistent
            end
            -model1d this-model-does-not-exist
            """);
        var code = new ScriptRunner(["OUT=results"], dryRun: false, keepGoing: false, CancellationToken.None).Run(Path.Combine(_dir, "run.tomo"));
        // The tolerated failure (a line starting with '-') is counted but does not stop the script.
        Assert.Equal(1, code);
        Assert.True(File.Exists(Path.Combine(_dir, "results", "ak135_d42.5.txt")));
        Assert.True(File.Exists(Path.Combine(_dir, "results", "iasp91_d42.5.txt")));
        Assert.True(File.Exists(Path.Combine(_dir, "results", "prem.txt")));
    }

    [Fact]
    public void AFailedCommandStopsTheScript()
    {
        File.WriteAllText(Path.Combine(_dir, "stop.tomo"), "model1d nothing-here\nmodel1d ak135 --out after.txt\n");
        Assert.Throws<ScriptFailedException>(() => new ScriptRunner([], false, false, CancellationToken.None).Run(Path.Combine(_dir, "stop.tomo")));
        Assert.False(File.Exists(Path.Combine(_dir, "after.txt")));
    }
}
