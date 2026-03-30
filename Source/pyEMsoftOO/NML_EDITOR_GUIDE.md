# EMsoftOO Namelist Editor — User Guide

A lightweight cross-platform GUI for browsing, editing, saving, and
running EMsoftOO namelist files.

---

## Requirements

- Python 3.9 or later
- tkinter (included with Python on macOS and Windows; on Ubuntu:
  `sudo apt install python3-tk`)
- The `emsoft` package installed (`pip install -e .` from the
  `Source/pyEMsoftOO` directory)
- A valid `~/.config/EMsoft/EMsoftConfig.json` configuration file
  (created by `EMsoftinit` or manually)

---

## Launching the Editor

The simplest way to launch the editor:

```bash
python3 -m emsoft.nml_editor
```

The editor reads `~/.config/EMsoft/EMsoftConfig.json` to find the
EMsoftOO installation path and locates the `NamelistTemplates/`
directory automatically. No command-line arguments are needed.

If the configuration file is not available, you can pass the template
directory explicitly:

```bash
python3 -m emsoft.nml_editor /path/to/EMsoftOO/NamelistTemplates
```

If neither method finds the templates, a folder chooser dialog will
appear on startup.

---

## Window Layout

The editor window has four main areas:

```
+---------------------------------------------------------------+
|  [Open .nml] [Save .nml] | [Set Work Dir] | [Run] [Stop]  Font|
+---------------+-------------------------------------------+---+
|               |                                            |   |
|  Template     |  Editor                                    |   |
|  List         |  (namelist file with syntax highlighting)  |   |
|               |                                            |   |
|  [Filter]     |                                            |   |
|               |                                            |   |
+---------------+-------------------------------------------+---+
|               |  Program Output                            |   |
|               |  (dark terminal-style panel)               |   |
|               |                        [Save .log] [Clear] |   |
+---------------+-------------------------------------------+---+
|  Work dir: /path/to/working/directory                         |
|  Status: Ready                                                |
+---------------------------------------------------------------+
```

### Template List (left panel)

- Shows all available namelist templates sorted alphabetically
- Click a template name to load it into the editor
- Use the **Filter** field at the top to narrow the list by typing
  part of a template name (e.g., typing "EBSD" shows only
  EBSD-related templates)

### Editor (upper right)

- Displays the namelist file content with syntax highlighting:
  - **Blue bold**: namelist group names (`&GroupName`) and terminators (`/`)
  - **Green**: comments (`! ...`)
  - **Dark blue**: parameter names (before `=`)
  - **Dark red**: string values (`'...'`)
  - **Purple**: boolean values (`.TRUE.`, `.FALSE.`)
  - **Bold green**: section dividers (`!====...`)
- Line numbers are shown in the left gutter
- Standard text editing: cut, copy, paste, undo/redo all work
- The editor tracks unsaved changes and warns before discarding them

### Program Output (lower right)

- Displays stdout and stderr from program execution in real time
- Dark terminal theme with color-coded output:
  - White: standard output
  - Red: error output (stderr)
  - Blue: status messages (command, working directory)
  - Green: completion messages

### Status Bar (bottom)

- Shows the current working directory
- Displays file status, save confirmations, and program run status

---

## Toolbar Buttons

### Open .nml

Opens a file dialog to load an existing `.nml` file from any location.
The dialog starts in the current working directory.

**Keyboard shortcut**: Cmd+O (macOS) / Ctrl+O (Windows/Linux)

### Save .nml

Saves the editor content as a `.nml` file. A file dialog appears with
the filename pre-filled (the template name with `.nml` extension). The
dialog defaults to the current working directory.

**Keyboard shortcut**: Cmd+S (macOS) / Ctrl+S (Windows/Linux)

### Set Work Dir

Opens a folder chooser to set the working directory. This directory is
used as:

- The default location for saving `.nml` files
- The working directory when running programs (the program's current
  directory will be set to this path)
- The default location for saving `.log` files
- The starting directory for the Open dialog

The current working directory is always shown in the status bar at the
bottom of the window.

### Run Program

Launches the EMsoftOO program corresponding to the current `.nml` file.
The program name is derived from the filename (e.g., `EMEBSDmaster.nml`
runs the `EMEBSDmaster` executable).

The editor finds executables by searching:

1. The system PATH
2. The `EMSOFTOO_BIN` environment variable
3. The `EMsoftLibraryLocation` path from `EMsoftConfig.json`

Before running:

- If the file has not been saved as `.nml` yet, you will be prompted to
  save it first
- If there are unsaved changes, you will be prompted to save them

**Keyboard shortcut**: Cmd+R (macOS) / Ctrl+R (Windows/Linux)

**Shift+Click**: Instead of running with the current editor file, a file
dialog opens to select any `.nml` file. The program name is still derived
from the selected filename. This is useful for running a previously saved
namelist without loading it into the editor. (Hover over the Run button
to see this hint.)

### Stop

Terminates a running program. Only active while a program is running.

### Font Size (- / +)

Adjusts the font size for the template list, editor, and output panel.
Range: 8 to 28 points.

**Keyboard shortcut**: Cmd+Plus / Cmd+Minus (macOS) / Ctrl+Plus /
Ctrl+Minus (Windows/Linux)

---

## Output Panel Buttons

### Save .log

Saves the program output to a `.log` file. The suggested filename
matches the `.nml` filename with a `.log` extension.

### Clear

Clears the output panel.

---

## Typical Workflow

### 1. Select a template

Click a template name in the left panel, e.g., **EMMCOpenCL**.

### 2. Edit the parameters

Modify the parameter values in the editor. The inline comments
describe each parameter. For example:

```
 &MCCLdata
 mode = 'full'
 xtalname = 'Ni.xtal',
 numsx = 501,
 sig = 70.0,
 ...
 /
```

### 3. Set the working directory

Click **Set Work Dir** and choose the folder where you want to run
the simulation (where input files like `.xtal` are located and where
output files will be written).

### 4. Save the namelist

Click **Save .nml** (or press Cmd/Ctrl+S). The file dialog defaults
to the working directory with the correct `.nml` filename.

### 5. Run the program

Click **Run Program** (or press Cmd/Ctrl+R). The program output
appears in the output panel in real time. When the program finishes,
a success or error message is shown.

### 6. Save the log (optional)

Click **Save .log** in the output panel to keep a record of the
program output.

### 7. Repeat

Select another template or modify the current namelist for another
run. The editor warns if you have unsaved changes before switching.

---

## Configuration File

The editor reads `~/.config/EMsoft/EMsoftConfig.json` which is the
standard EMsoftOO configuration file. The relevant fields are:

| Field | Used For |
|-------|----------|
| `EMsoftpathname` | Locating the `NamelistTemplates/` directory |
| `EMsoftLibraryLocation` | Finding EMsoftOO executables for the Run button |

Example `EMsoftConfig.json`:

```json
{
    "EMsoftpathname": "/Users/mdg/Files/EMsoftOO/",
    "EMsoftLibraryLocation": "/Users/mdg/Files/EMsoftOOBuild/Bin/",
    "EMdatapathname": "/Users/mdg/Files/EMPlay/",
    ...
}
```

This file is created by the `EMsoftinit` program or can be edited
manually.

---

## Platform Notes

### macOS

- Keyboard shortcuts use Cmd (Command) key
- tkinter buttons may occasionally require a precise click; this is a
  known macOS tkinter limitation
- The editor uses the Menlo font for the editor and output panels

### Windows

- Keyboard shortcuts use Ctrl key
- High-DPI displays are supported automatically
- The editor uses the Consolas font

### Linux (Ubuntu)

- Keyboard shortcuts use Ctrl key
- If tkinter is not installed: `sudo apt install python3-tk`
- The editor uses the system Monospace font

---

## Troubleshooting

### "Could not find NamelistTemplates directory"

Make sure `~/.config/EMsoft/EMsoftConfig.json` exists and contains a
valid `EMsoftpathname` entry pointing to the EMsoftOO source directory.
Alternatively, pass the path as a command-line argument:

```bash
python3 -m emsoft.nml_editor /path/to/EMsoftOO/NamelistTemplates
```

### Run button is grayed out

The Run button is enabled only when:

1. The current file is saved as a `.nml` file (not a `.template`)
2. The corresponding executable is found (via PATH, EMSOFTOO_BIN,
   or EMsoftConfig.json)

Check that the program name matches the filename (e.g., `EMEBSDmaster.nml`
requires the `EMEBSDmaster` executable) and that the executable directory
is correctly configured.

### Program does not produce output

Some EMsoftOO programs use interactive terminal prompts (e.g., `EMmkxtal`).
These programs are not compatible with the Run button since they require
keyboard input during execution. The Run feature works best with programs
that read all parameters from the `.nml` file.
