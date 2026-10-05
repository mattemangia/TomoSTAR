// Copyright 2026 Matteo Mangiagalli
// SPDX-License-Identifier: Apache-2.0

using TomoStar.Core.IO;

namespace TomoStar.Cli;

/// <summary>A command: what it does, how it is called, and the code that runs it.</summary>
public sealed record CommandInfo(string Name, string Summary, string Usage, Func<CommandContext, int> Run);

/// <summary>
/// The entry point of <c>tomostar</c>: finds the command, builds its context (arguments,
/// configuration, data) and runs it. Scripts call <see cref="Execute"/> for each of their lines, so a
/// command behaves the same typed at the prompt or written in a script.
/// </summary>
public static class Program
{
    public static readonly IReadOnlyList<CommandInfo> Commands =
    [
        new("synth", "Make a synthetic data set with a known answer (stations, events, picks, t*, waveforms).",
            "synth --out DIR [--with-waveforms] [--events N] [--stations N] [--seed N] [--set Synthetic.X=...]", SetupCommands.Synth),
        new("grid", "Propose the inversion grid and the 1-D profile on it from the data (or from --bounds).",
            "grid DATA --out DIR [--h KM] [--dz KM] [--margin KM] [--top KM] [--bottom KM] [--bounds LON1,LON2,LAT1,LAT2] [--model M]", SetupCommands.Grid),
        new("model1d", "List the published 1-D models, show or export one, choose one for an area, sample it on a grid.",
            "model1d --list | model1d NAME|FILE [--out FILE] [--grid grid.json [--profile FILE]] | model1d --area LON1,LON2,LAT1,LAT2", SetupCommands.Model1D),
        new("pick", "Automatic P and S picks (STA/LTA trigger, AIC onset) around the arrivals the model predicts.",
            "pick DATA --waveforms DIR --out DIR [--only-automatic] [--vp V.qvol --vs V.qvol]", PipelineCommands.Pick),
        new("relocate", "Relocate the events in a 1-D or 3-D model: absolute (grid search + Geiger) or double difference.",
            "relocate DATA --out DIR [--method absolute|dd] [--fix-depth] [--vp V.qvol --vs V.qvol]", PipelineCommands.Relocate),
        new("minimum1d", "Minimum 1-D model: layered inversion with hypocentres and station terms (Kissling et al. 1994).",
            "minimum1d DATA --out DIR [--iterations N] [--damping D] [--smoothing S]", InversionCommands.Minimum1D),
        new("invert", "Travel-time tomography: Vp and Vs (or Vp/Vs) jointly with hypocentres and station terms.",
            "invert DATA --out DIR [--iterations N] [--damping D] [--smoothing S] [--adaptive] [--layered] [--vp V.qvol --vs V.qvol]", InversionCommands.Invert),
        new("lcurve", "Trade-off (L-curve) of damping and smoothing for the velocity or the Q tomography; recommends a pair.",
            "lcurve DATA --out DIR [--kind velocity|q] [--dampings 1,2,5,...] [--smoothings 5,20,80] [--save-models] [--update-config FILE]", InversionCommands.LCurve),
        new("resolution", "Checkerboard, spike and tabular-body tests for Vp, Vs, Vp/Vs, their leakage, and Q.",
            "resolution DATA --out DIR [--pattern checkerboard|spike|body] [--target Vp|Vs|VpVs|Leakage|Q|QLeakage] [--cell KM] [--cell-depth KM] [--amplitude %] [--spike LON,LAT,Z,R,%] [--body LON,LAT,Z,STRIKE,DIP,L,W,T,%]", InversionCommands.Resolution),
        new("tstar", "Measure t* from the displacement spectra of the P (and S) waves.",
            "tstar DATA --waveforms DIR --out DIR [--stationxml FILE]", PipelineCommands.TStar),
        new("qtomo", "Attenuation (Qp or Qs) tomography from t*, along the rays of the velocity model.",
            "qtomo DATA --out DIR --vp Vp.qvol --vs Vs.qvol [--damping D] [--smoothing S] [--phase P|S]", InversionCommands.QTomo),
        new("export", "Convert .qvol volumes to CSV or VTK (the layouts QUIVER imports).",
            "export FILE.qvol [...] [--csv] [--vtk] [--out DIR]", SetupCommands.Export),
        new("register", "Add a run folder and its volumes to a QUIVER project.",
            "register RUN_FOLDER --project QUIVER_PROJECT [--name NAME]", SetupCommands.Register),
        new("config", "Show the effective configuration, or write a template with every setting (--template).",
            "config [--template --out FILE] | config DATA [--config FILE] [--set Path=value]", SetupCommands.Config),
        new("info", "Version, CPU and SIMD width, OpenCL devices and their self-tests.", "info", SetupCommands.Info),
    ];

    private const string GlobalOptions = """
        Data (any command that reads data):
          --data PATH           folder with stations.csv, events.csv, picks.csv, tstar.csv, or a QUIVER project
          --stations F --events F --picks F --tstar F   single tables (replace those of --data)
          --quakeml F --stationxml F                    add events and picks, or stations and responses
          --hypocentres F       event table whose locations replace those of the same events
          --catalog-locations   QUIVER project: start from the catalogue locations, not the relocations
        Model and grid:
          --grid F              grid.json (from 'tomostar grid') or a QUIVER project
          --model M             1-D model: library name, file of 'depth vp vs' lines, auto
          --vp F --vs F         3-D velocities (.qvol) instead of the 1-D model
        Settings and output:
          --config F            JSON configuration ('tomostar config --template' writes one)
          --set Path=value      override one setting, e.g. --set Tomography.Smoothing=30 (repeatable)
          --out DIR  --name N   run folder and run name (default runs/<date_time>_<command>)
          --formats csv,vtk     volume formats besides .qvol (also --csv, --vtk)
          --register PROJECT    copy the results into a QUIVER project ('data' = the project given as --data)
          --adaptive            adaptive (octree) parameterisation; --straight: straight rays
          --no-opencl --threads N --quiet --keep-work
        """;

    public static int Main(string[] args)
    {
        using var cts = new CancellationTokenSource();
        Console.CancelKeyPress += (_, e) =>
        {
            e.Cancel = true;
            cts.Cancel();
            Console.Error.WriteLine("Stopping at the next checkpoint...");
        };
        return Execute(args, cts.Token);
    }

    /// <summary>Runs one command line (without the program name); returns the exit code.</summary>
    public static int Execute(IReadOnlyList<string> args, CancellationToken ct)
    {
        if (args.Count == 0 || args[0] is "help" or "--help" or "-h")
        {
            Help(args.Count > 1 ? args[1] : null);
            return 0;
        }
        if (args[0] is "version" or "--version")
        {
            Console.WriteLine($"TomoSTAR {typeof(Program).Assembly.GetName().Version?.ToString(3)}");
            return 0;
        }
        CommandContext? context = null;
        try
        {
            if (args[0] == "run") return RunScript(args.Skip(1).ToList(), ct);
            var command = Commands.FirstOrDefault(c => c.Name.Equals(args[0], StringComparison.OrdinalIgnoreCase))
                          ?? throw new UsageException($"Unknown command '{args[0]}'.");
            var parsed = Arguments.Parse(args.Skip(1).ToList());
            if (parsed.Flag("help"))
            {
                Help(command.Name);
                return 0;
            }
            // A positional argument names the data for the commands that read data.
            if (command.Name is not ("model1d" or "export" or "register") && parsed.Positional.Count > 0 && parsed.Get("data") == null)
                parsed = Arguments.Parse(args.Skip(1).Where(a => a != parsed.Positional[0]).Concat(["--data", parsed.Positional[0]]).ToList());
            context = new CommandContext(command.Name, parsed, ct);
            var code = command.Run(context);
            foreach (var u in parsed.Unused()) Console.Error.WriteLine($"Warning: option {u} was not used by '{command.Name}'.");
            return code;
        }
        catch (UsageException ex)
        {
            Console.Error.WriteLine($"Error: {ex.Message}");
            Console.Error.WriteLine($"See 'tomostar help{(args.Count > 0 && args[0] != "run" ? " " + args[0] : "")}'.");
            return 2;
        }
        catch (ScriptFailedException ex)
        {
            Console.Error.WriteLine(ex.Message);
            return 1;
        }
        catch (OperationCanceledException)
        {
            Console.Error.WriteLine("Cancelled.");
            return 130;
        }
        catch (Exception ex) when (ex is IOException or InvalidDataException or InvalidOperationException or ArgumentException
                                       or KeyNotFoundException or UnauthorizedAccessException or System.Text.Json.JsonException or FormatException)
        {
            Console.Error.WriteLine($"Error: {ex.Message}");
            if (Environment.GetEnvironmentVariable("TOMOSTAR_DEBUG") == "1") Console.Error.WriteLine(ex);
            if (context != null)
            {
                try { RunWriter.WriteLog(context.RunFolder, [.. context.LogLines, "Error: " + ex.Message]); } catch (IOException) { }
            }
            return 1;
        }
        finally
        {
            context?.CleanUp();
        }
    }

    private static int RunScript(List<string> args, CancellationToken ct)
    {
        var parsed = Arguments.Parse(args);
        if (parsed.Positional.Count != 1) throw new UsageException("run SCRIPT [--var NAME=value ...] [--dry-run] [--keep-going]");
        var runner = new ScriptRunner(parsed.All("var"), parsed.Flag("dry-run"), parsed.Flag("keep-going"), ct);
        return runner.Run(parsed.Positional[0]);
    }

    private static void Help(string? command)
    {
        if (command != null && Commands.FirstOrDefault(c => c.Name.Equals(command, StringComparison.OrdinalIgnoreCase)) is { } info)
        {
            Console.WriteLine($"tomostar {info.Name}: {info.Summary}\n\nUsage: tomostar {info.Usage}\n");
            Console.WriteLine(GlobalOptions);
            return;
        }
        if (command == "run")
        {
            Console.WriteLine("""
                tomostar run SCRIPT [--var NAME=value ...] [--dry-run] [--keep-going]

                Runs a TomoSTAR script: one command per line (without 'tomostar'), plus
                  set NAME = value | default NAME = value | config FILE | options --opt ... | echo TEXT
                  include FILE | foreach NAME in A B C ... end | if exists PATH ... end | if missing PATH ... end | exit
                ${NAME}, ${env:NAME} and ${json:FILE:Path.To.Value} are expanded; a line starting with '-' may fail
                without stopping the script; relative paths are relative to the script's folder.
                See examples/*.tomo and the README.
                """);
            return;
        }
        Console.WriteLine("""
            TomoSTAR: seismic travel-time and attenuation tomography from the command line.

            Usage: tomostar COMMAND [DATA] [OPTIONS]    tomostar help COMMAND    tomostar run SCRIPT
            """);
        Console.WriteLine("Commands:");
        foreach (var c in Commands) Console.WriteLine($"  {c.Name,-11} {c.Summary}");
        Console.WriteLine($"  {"run",-11} Run a script of commands (tomostar help run).");
        Console.WriteLine();
        Console.WriteLine(GlobalOptions);
    }
}
