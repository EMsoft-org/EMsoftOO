"""
EMsoftOO Namelist Editor — a lightweight GUI for editing namelist template files.

Launch with:
    python -m emsoft.nml_editor
    python -m emsoft.nml_editor /path/to/NamelistTemplates
"""

import os
import sys
import platform
import shutil
import subprocess
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

# Detect platform for keyboard shortcuts
IS_MAC = platform.system() == 'Darwin'
MOD_KEY = 'Command' if IS_MAC else 'Control'
MOD_LABEL = 'Cmd' if IS_MAC else 'Ctrl'


def find_templates_dir(explicit_path=None):
    """Locate the NamelistTemplates directory."""
    candidates = []

    if explicit_path:
        candidates.append(explicit_path)

    env = os.environ.get('EMSOFTOO_TEMPLATES')
    if env:
        candidates.append(env)

    # Walk up from this file to find the repo root
    here = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        candidate = os.path.join(here, 'NamelistTemplates')
        candidates.append(candidate)
        here = os.path.dirname(here)

    for path in candidates:
        if os.path.isdir(path):
            templates = [f for f in os.listdir(path) if f.endswith('.template')]
            if templates:
                return path

    return None


def _find_executable(program_name):
    """Find an EMsoftOO executable by name."""
    # 1. Check PATH
    path = shutil.which(program_name)
    if path:
        return path

    # 2. Check EMSOFTOO_BIN environment variable
    bin_dir = os.environ.get('EMSOFTOO_BIN')
    if bin_dir:
        candidate = os.path.join(bin_dir, program_name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate

    return None


def _program_name_from_file(filepath):
    """Derive the EMsoftOO program name from a .nml or .template filename."""
    base = os.path.basename(filepath)
    name = base.replace('.nml', '').replace('.template', '')
    # Template names like "BetheParameters" don't correspond to programs
    # Program names start with "EM" (convention)
    if name.startswith('EM'):
        return name
    return None


class NmlEditor:
    """Main application window."""

    def __init__(self, root, templates_dir):
        self.root = root
        self.templates_dir = templates_dir
        self.current_file = None
        self.modified = False
        self._highlight_job = None
        self._process = None
        self.work_dir = os.getcwd()

        self.root.title('EMsoftOO Namelist Editor')
        self.root.geometry('1000x750')
        self.root.minsize(700, 500)

        # Handle window close
        self.root.protocol('WM_DELETE_WINDOW', self.on_close)

        # Configure editor font (platform-specific)
        if IS_MAC:
            self.editor_font = ('Menlo', 13)
        elif platform.system() == 'Windows':
            self.editor_font = ('Consolas', 11)
        else:
            self.editor_font = ('Monospace', 11)

        self._build_ui()
        self._bind_shortcuts()
        self._apply_theme()

    def _build_ui(self):
        """Create all UI elements."""
        # --- Top toolbar ---
        toolbar = tk.Frame(self.root, padx=5, pady=5)
        toolbar.pack(fill=tk.X)

        self.open_btn = tk.Button(toolbar, text='  Open .nml  ', command=self.open_file)
        self.open_btn.pack(side=tk.LEFT, padx=4)

        self.save_btn = tk.Button(toolbar, text=f'  Save .nml ({MOD_LABEL}+S)  ',
                                  command=self.save_file)
        self.save_btn.pack(side=tk.LEFT, padx=4)

        tk.Label(toolbar, text='  |  ').pack(side=tk.LEFT)

        self.workdir_btn = tk.Button(toolbar, text='  Set Work Dir  ',
                                     command=self.set_work_dir)
        self.workdir_btn.pack(side=tk.LEFT, padx=4)

        tk.Label(toolbar, text='  |  ').pack(side=tk.LEFT)

        self.run_btn = tk.Button(toolbar, text='  Run Program  ', command=self.run_program,
                                 state=tk.DISABLED)
        self.run_btn.pack(side=tk.LEFT, padx=4)

        self.stop_btn = tk.Button(toolbar, text='  Stop  ', command=self.stop_program,
                                  state=tk.DISABLED, fg='red')
        self.stop_btn.pack(side=tk.LEFT, padx=4)

        # Font size controls (right side of toolbar)
        font_frame = tk.Frame(toolbar)
        font_frame.pack(side=tk.RIGHT, padx=4)
        tk.Label(font_frame, text='Font:').pack(side=tk.LEFT)
        tk.Button(font_frame, text=' \u2212 ', command=self._font_smaller).pack(side=tk.LEFT, padx=1)
        self.font_size_var = tk.StringVar(value=str(self.editor_font[1]))
        tk.Label(font_frame, textvariable=self.font_size_var, width=3,
                 anchor=tk.CENTER).pack(side=tk.LEFT)
        tk.Button(font_frame, text=' + ', command=self._font_larger).pack(side=tk.LEFT, padx=1)

        # --- Main paned layout: template list on left, editor+output on right ---
        h_paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        h_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=(0, 0))

        # --- Left panel: template list with search ---
        left_frame = ttk.Frame(h_paned, padding=2)
        h_paned.add(left_frame, weight=0)

        ttk.Label(left_frame, text='Templates:').pack(anchor=tk.W)

        # Search entry
        search_frame = ttk.Frame(left_frame)
        search_frame.pack(fill=tk.X, pady=(2, 4))
        ttk.Label(search_frame, text='Filter:').pack(side=tk.LEFT)
        self.search_var = tk.StringVar()
        self.search_var.trace_add('write', self._on_search_changed)
        search_entry = ttk.Entry(search_frame, textvariable=self.search_var, width=20)
        search_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 0))

        # Template listbox with scrollbar
        list_frame = ttk.Frame(left_frame)
        list_frame.pack(fill=tk.BOTH, expand=True)

        list_scroll = ttk.Scrollbar(list_frame, orient=tk.VERTICAL)
        list_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        self.template_list = tk.Listbox(list_frame, width=28, activestyle='dotbox',
                                         yscrollcommand=list_scroll.set,
                                         font=('TkDefaultFont', 11))
        self.template_list.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        list_scroll.config(command=self.template_list.yview)

        self.template_list.bind('<<ListboxSelect>>', self.on_template_selected)

        # Populate the list
        self.all_templates = sorted(
            f for f in os.listdir(self.templates_dir) if f.endswith('.template')
        )
        self.display_names = [f.replace('.template', '') for f in self.all_templates]
        self._populate_list(self.display_names)

        # --- Right panel: editor (top) + output (bottom) ---
        right_frame = ttk.Frame(h_paned, padding=2)
        h_paned.add(right_frame, weight=1)

        v_paned = ttk.PanedWindow(right_frame, orient=tk.VERTICAL)
        v_paned.pack(fill=tk.BOTH, expand=True)

        # --- Editor pane ---
        editor_container = ttk.Frame(v_paned)
        v_paned.add(editor_container, weight=3)

        editor_frame = ttk.Frame(editor_container)
        editor_frame.pack(fill=tk.BOTH, expand=True)

        # Line numbers
        self.linenums = tk.Text(editor_frame, width=4, padx=4, pady=4,
                                takefocus=0, border=0, state='disabled',
                                background='#f0f0f0', foreground='#999999',
                                font=self.editor_font)
        self.linenums.pack(side=tk.LEFT, fill=tk.Y)

        # Editor scrollbar
        scrollbar = ttk.Scrollbar(editor_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Main text editor
        self.editor = tk.Text(editor_frame, wrap=tk.NONE, undo=True,
                              padx=6, pady=4, insertwidth=2,
                              font=self.editor_font,
                              yscrollcommand=self._on_scroll)
        self.editor.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self._scroll_both)

        # Track modifications
        self.editor.bind('<<Modified>>', self._on_modified)
        self.editor.bind('<KeyRelease>', self._schedule_highlight)

        # --- Output pane ---
        output_container = ttk.Frame(v_paned)
        v_paned.add(output_container, weight=2)

        # Output header with label and Save .log button
        output_header = tk.Frame(output_container, padx=2, pady=2)
        output_header.pack(fill=tk.X)
        tk.Label(output_header, text='Program Output:',
                 font=('TkDefaultFont', 10, 'bold')).pack(side=tk.LEFT)
        tk.Button(output_header, text=' Clear ', command=self.clear_output).pack(side=tk.RIGHT, padx=2)
        tk.Button(output_header, text=' Save .log ', command=self.save_log).pack(side=tk.RIGHT, padx=2)

        # Output text widget
        output_frame = ttk.Frame(output_container)
        output_frame.pack(fill=tk.BOTH, expand=True)

        output_scroll = ttk.Scrollbar(output_frame)
        output_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        self.output = tk.Text(output_frame, wrap=tk.WORD, state=tk.DISABLED,
                              padx=6, pady=4, font=self.editor_font,
                              background='#1e1e1e', foreground='#cccccc',
                              insertbackground='#cccccc',
                              yscrollcommand=output_scroll.set)
        self.output.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        output_scroll.config(command=self.output.yview)

        # Output text tags
        self.output.tag_configure('stderr', foreground='#f44747')
        self.output.tag_configure('info', foreground='#569cd6')
        self.output.tag_configure('success', foreground='#4ec9b0')

        # --- Status bar (two rows: work dir + status) ---
        status_frame = tk.Frame(self.root)
        status_frame.pack(fill=tk.X, side=tk.BOTTOM)

        self.workdir_var = tk.StringVar(value=f'Work dir: {self.work_dir}')
        workdir_label = ttk.Label(status_frame, textvariable=self.workdir_var,
                                  relief=tk.SUNKEN, anchor=tk.W, padding=(5, 2),
                                  foreground='#555555')
        workdir_label.pack(fill=tk.X)

        self.status_var = tk.StringVar(value='Select a template from the list to begin')
        status = ttk.Label(status_frame, textvariable=self.status_var,
                           relief=tk.SUNKEN, anchor=tk.W, padding=(5, 2))
        status.pack(fill=tk.X)

        # --- Syntax highlighting tags ---
        self.editor.tag_configure('comment', foreground='#6a9955')
        self.editor.tag_configure('section', foreground='#2e8b57',
                                  font=self.editor_font + ('bold',))
        self.editor.tag_configure('group', foreground='#0000cc',
                                  font=self.editor_font + ('bold',))
        self.editor.tag_configure('param', foreground='#001080')
        self.editor.tag_configure('string', foreground='#a31515')
        self.editor.tag_configure('boolean', foreground='#7b2fa0')

    def _populate_list(self, names):
        """Fill the listbox with template names."""
        self.template_list.delete(0, tk.END)
        for name in names:
            self.template_list.insert(tk.END, name)

    def _on_search_changed(self, *args):
        """Filter the template list based on search text."""
        query = self.search_var.get().lower()
        if query:
            filtered = [n for n in self.display_names if query in n.lower()]
        else:
            filtered = self.display_names
        self._populate_list(filtered)

    def _font_smaller(self):
        """Decrease font size by 1 (minimum 8)."""
        size = self.editor_font[1]
        if size > 8:
            self._set_font_size(size - 1)

    def _font_larger(self):
        """Increase font size by 1 (maximum 28)."""
        size = self.editor_font[1]
        if size < 28:
            self._set_font_size(size + 1)

    def _set_font_size(self, size):
        """Update font size across all widgets."""
        self.editor_font = (self.editor_font[0], size)
        self.font_size_var.set(str(size))

        self.editor.configure(font=self.editor_font)
        self.linenums.configure(font=self.editor_font)
        self.output.configure(font=self.editor_font)
        self.template_list.configure(font=('TkDefaultFont', size))

        # Re-apply highlighting tags with new font size
        self.editor.tag_configure('section',
                                  font=self.editor_font + ('bold',))
        self.editor.tag_configure('group',
                                  font=self.editor_font + ('bold',))

        self._update_line_numbers()

    def _bind_shortcuts(self):
        """Set up keyboard shortcuts."""
        mod = 'Command' if IS_MAC else 'Control'
        self.root.bind(f'<{mod}-s>', lambda e: self.save_file())
        self.root.bind(f'<{mod}-S>', lambda e: self.save_file())
        self.root.bind(f'<{mod}-o>', lambda e: self.open_file())
        self.root.bind(f'<{mod}-O>', lambda e: self.open_file())
        self.root.bind(f'<{mod}-equal>', lambda e: self._font_larger())
        self.root.bind(f'<{mod}-plus>', lambda e: self._font_larger())
        self.root.bind(f'<{mod}-minus>', lambda e: self._font_smaller())
        self.root.bind(f'<{mod}-r>', lambda e: self.run_program())
        self.root.bind(f'<{mod}-R>', lambda e: self.run_program())

    def _apply_theme(self):
        """Apply consistent styling."""
        style = ttk.Style()
        try:
            if IS_MAC:
                style.theme_use('aqua')
            elif platform.system() == 'Windows':
                style.theme_use('vista')
            else:
                style.theme_use('clam')
        except tk.TclError:
            pass

    # --- Scrolling ---

    def _scroll_both(self, *args):
        self.editor.yview(*args)
        self.linenums.yview(*args)

    def _on_scroll(self, first, last):
        self.linenums.yview_moveto(first)
        for child in self.editor.master.winfo_children():
            if isinstance(child, ttk.Scrollbar):
                child.set(first, last)
                break

    # --- Line numbers ---

    def _update_line_numbers(self):
        self.linenums.config(state='normal')
        self.linenums.delete('1.0', tk.END)

        line_count = int(self.editor.index('end-1c').split('.')[0])
        lines = '\n'.join(str(i) for i in range(1, line_count + 1))
        self.linenums.insert('1.0', lines)

        width = max(4, len(str(line_count)) + 1)
        self.linenums.config(width=width, state='disabled')

    # --- Syntax highlighting ---

    def _schedule_highlight(self, event=None):
        """Debounce highlighting to avoid lag on fast typing."""
        if self._highlight_job:
            self.root.after_cancel(self._highlight_job)
        self._highlight_job = self.root.after(100, self._apply_highlighting)

    def _apply_highlighting(self):
        """Apply syntax highlighting to the entire editor content."""
        self._highlight_job = None

        for tag in ('comment', 'section', 'group', 'param', 'string', 'boolean'):
            self.editor.tag_remove(tag, '1.0', tk.END)

        content = self.editor.get('1.0', tk.END)
        lines = content.split('\n')

        for i, line in enumerate(lines):
            line_start = f'{i + 1}.0'
            line_end = f'{i + 1}.end'
            stripped = line.lstrip()

            if not stripped:
                continue

            if stripped.startswith('&'):
                self.editor.tag_add('group', line_start, line_end)
                continue

            if stripped == '/':
                self.editor.tag_add('group', line_start, line_end)
                continue

            if stripped.startswith('!') and '===' in stripped:
                self.editor.tag_add('section', line_start, line_end)
                continue

            if stripped.startswith('!'):
                self.editor.tag_add('comment', line_start, line_end)
                continue

            if '=' in line and not stripped.startswith('!'):
                eq_pos = line.index('=')
                param_end = f'{i + 1}.{eq_pos}'
                self.editor.tag_add('param', line_start, param_end)

                value_part = line[eq_pos + 1:]
                in_string = False
                comment_offset = None
                for j, ch in enumerate(value_part):
                    if ch == "'" and not in_string:
                        in_string = True
                    elif ch == "'" and in_string:
                        in_string = False
                    elif ch == '!' and not in_string:
                        comment_offset = eq_pos + 1 + j
                        break

                if comment_offset is not None:
                    self.editor.tag_add('comment',
                                        f'{i + 1}.{comment_offset}', line_end)
                    value_text = line[eq_pos + 1:comment_offset]
                else:
                    value_text = value_part

                val_start = eq_pos + 1
                in_str = False
                str_begin = 0
                for j, ch in enumerate(value_text):
                    if ch == "'" and not in_str:
                        in_str = True
                        str_begin = val_start + j
                    elif ch == "'" and in_str:
                        in_str = False
                        self.editor.tag_add('string',
                                            f'{i + 1}.{str_begin}',
                                            f'{i + 1}.{val_start + j + 1}')

                val_upper = value_text.upper()
                for bval in ['.TRUE.', '.FALSE.']:
                    idx = val_upper.find(bval)
                    if idx >= 0:
                        self.editor.tag_add('boolean',
                                            f'{i + 1}.{val_start + idx}',
                                            f'{i + 1}.{val_start + idx + len(bval)}')

        self._update_line_numbers()

    # --- Modification tracking ---

    def _on_modified(self, event=None):
        if self.editor.edit_modified():
            if not self.modified:
                self.modified = True
                self._update_title()
            self.editor.edit_modified(False)

    def _update_title(self):
        title = 'EMsoftOO Namelist Editor'
        if self.current_file:
            name = os.path.basename(self.current_file)
            title += f' \u2014 {name}'
        if self.modified:
            title += ' (modified)'
        self.root.title(title)

    def _update_run_button(self):
        """Enable/disable the Run button based on current file."""
        if self.current_file and self.current_file.endswith('.nml'):
            prog = _program_name_from_file(self.current_file)
            if prog and _find_executable(prog):
                self.run_btn.config(state=tk.NORMAL)
                return
        self.run_btn.config(state=tk.DISABLED)

    # --- Template loading ---

    def on_template_selected(self, event=None):
        selection = self.template_list.curselection()
        if not selection:
            return

        name = self.template_list.get(selection[0])

        if self.modified:
            if not messagebox.askyesno('Unsaved Changes',
                    'Current file has unsaved changes. Discard them?'):
                return

        template_file = name + '.template'
        filepath = os.path.join(self.templates_dir, template_file)

        if os.path.isfile(filepath):
            self._load_file(filepath)
            line_count = int(self.editor.index('end-1c').split('.')[0])
            self.status_var.set(f'{template_file} loaded ({line_count} lines)')
        else:
            self.status_var.set(f'Template not found: {template_file}')

    def _load_file(self, filepath):
        """Load a file into the editor."""
        with open(filepath, 'r') as f:
            content = f.read()

        self.editor.delete('1.0', tk.END)
        self.editor.insert('1.0', content)

        self.editor.edit_modified(False)
        self.editor.edit_reset()
        self.modified = False
        self.current_file = filepath
        self._update_title()
        self._apply_highlighting()
        self._update_run_button()
        self.editor.see('1.0')

    # --- Working directory ---

    def set_work_dir(self):
        """Choose the working directory for saving .nml files and running programs."""
        d = filedialog.askdirectory(
            title='Select Working Directory',
            initialdir=self.work_dir)
        if d:
            self.work_dir = d
            self.workdir_var.set(f'Work dir: {self.work_dir}')
            self.status_var.set(f'Working directory set to: {d}')

    # --- File operations ---

    def open_file(self):
        """Open an existing .nml or .template file."""
        if self.modified:
            if not messagebox.askyesno('Unsaved Changes',
                    'Current file has unsaved changes. Discard them?'):
                return

        filepath = filedialog.askopenfilename(
            title='Open Namelist File',
            filetypes=[
                ('Namelist files', '*.nml'),
                ('Template files', '*.template'),
                ('All files', '*.*'),
            ],
            initialdir=self.work_dir,
        )
        if filepath:
            self._load_file(filepath)
            name = os.path.basename(filepath)
            line_count = int(self.editor.index('end-1c').split('.')[0])
            self.status_var.set(f'{name} loaded ({line_count} lines)')

    def save_file(self):
        """Save the editor content as a .nml file."""
        if self.current_file:
            base = os.path.basename(self.current_file)
            default_name = base.replace('.template', '.nml')
            if not default_name.endswith('.nml'):
                default_name = base
        else:
            default_name = 'unnamed.nml'

        filepath = filedialog.asksaveasfilename(
            title='Save Namelist File',
            defaultextension='.nml',
            filetypes=[
                ('Namelist files', '*.nml'),
                ('All files', '*.*'),
            ],
            initialfile=default_name,
            initialdir=self.work_dir,
        )
        if filepath:
            content = self.editor.get('1.0', 'end-1c')
            with open(filepath, 'w') as f:
                f.write(content)
                if not content.endswith('\n'):
                    f.write('\n')

            self.current_file = filepath
            self.modified = False
            self._update_title()
            self._update_run_button()
            self.status_var.set(f'Saved: {filepath}')

    # --- Program execution ---

    def run_program(self):
        """Launch the EMsoftOO program for the current .nml file."""
        if self._process is not None:
            messagebox.showinfo('Already Running', 'A program is already running.')
            return

        if not self.current_file or not self.current_file.endswith('.nml'):
            # Offer to save first
            if messagebox.askyesno('Save First',
                    'The file must be saved as .nml before running.\nSave now?'):
                self.save_file()
                if not self.current_file or not self.current_file.endswith('.nml'):
                    return
            else:
                return

        if self.modified:
            if messagebox.askyesno('Unsaved Changes',
                    'Save changes before running?'):
                # Quick-save to current path
                content = self.editor.get('1.0', 'end-1c')
                with open(self.current_file, 'w') as f:
                    f.write(content)
                    if not content.endswith('\n'):
                        f.write('\n')
                self.modified = False
                self._update_title()

        prog_name = _program_name_from_file(self.current_file)
        if not prog_name:
            messagebox.showerror('Error',
                'Could not determine program name from filename.\n'
                'EMsoftOO programs start with "EM".')
            return

        exe_path = _find_executable(prog_name)
        if not exe_path:
            messagebox.showerror('Program Not Found',
                f'Could not find executable: {prog_name}\n\n'
                f'Make sure it is in your PATH or set the\n'
                f'EMSOFTOO_BIN environment variable to the\n'
                f'directory containing EMsoftOO executables.')
            return

        nml_path = os.path.abspath(self.current_file)
        work_dir = self.work_dir
        nml_basename = os.path.basename(nml_path)

        # Clear output and show start message
        self.clear_output()
        self._append_output(f'Running: {prog_name} {nml_basename}\n', 'info')
        self._append_output(f'Working directory: {work_dir}\n', 'info')
        self._append_output(f'{"-" * 60}\n', 'info')

        self.status_var.set(f'Running {prog_name}...')
        self.run_btn.config(state=tk.DISABLED)
        self.stop_btn.config(state=tk.NORMAL)

        # Launch the process in a background thread
        self._process_thread = threading.Thread(
            target=self._run_in_thread,
            args=(exe_path, nml_basename, work_dir),
            daemon=True)
        self._process_thread.start()

    def _run_in_thread(self, exe_path, nml_file, work_dir):
        """Run the program in a background thread, streaming output to the GUI."""
        try:
            self._process = subprocess.Popen(
                [exe_path, nml_file],
                cwd=work_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1)

            # Read stdout and stderr in parallel
            def read_stream(stream, tag):
                try:
                    for line in stream:
                        self.root.after(0, self._append_output, line, tag)
                except Exception:
                    pass

            stderr_thread = threading.Thread(
                target=read_stream, args=(self._process.stderr, 'stderr'),
                daemon=True)
            stderr_thread.start()

            # Read stdout in this thread
            read_stream(self._process.stdout, None)

            # Wait for process to finish
            self._process.wait()
            stderr_thread.join(timeout=2)

            retcode = self._process.returncode
            self.root.after(0, self._on_process_finished, retcode)

        except Exception as e:
            self.root.after(0, self._append_output,
                           f'\nError launching program: {e}\n', 'stderr')
            self.root.after(0, self._on_process_finished, -1)

    def _on_process_finished(self, retcode):
        """Called when the program finishes (on the main thread)."""
        self._process = None
        self.run_btn.config(state=tk.NORMAL)
        self.stop_btn.config(state=tk.DISABLED)
        self._update_run_button()

        self._append_output(f'\n{"-" * 60}\n', 'info')
        if retcode == 0:
            self._append_output('Program completed successfully.\n', 'success')
            self.status_var.set('Program completed successfully')
        else:
            self._append_output(f'Program exited with code {retcode}.\n', 'stderr')
            self.status_var.set(f'Program exited with code {retcode}')

    def stop_program(self):
        """Stop the running program."""
        if self._process is not None:
            self._process.terminate()
            self._append_output('\n*** Program terminated by user ***\n', 'stderr')
            self.status_var.set('Program terminated')

    # --- Output pane ---

    def _append_output(self, text, tag=None):
        """Append text to the output pane (must be called from main thread)."""
        self.output.config(state=tk.NORMAL)
        if tag:
            self.output.insert(tk.END, text, tag)
        else:
            self.output.insert(tk.END, text)
        self.output.see(tk.END)
        self.output.config(state=tk.DISABLED)

    def clear_output(self):
        """Clear the output pane."""
        self.output.config(state=tk.NORMAL)
        self.output.delete('1.0', tk.END)
        self.output.config(state=tk.DISABLED)

    def save_log(self):
        """Save the output pane content as a .log file."""
        content = self.output.get('1.0', 'end-1c')
        if not content.strip():
            messagebox.showinfo('Empty Output', 'There is no output to save.')
            return

        if self.current_file:
            default_name = os.path.basename(self.current_file).replace('.nml', '.log')
        else:
            default_name = 'output.log'

        filepath = filedialog.asksaveasfilename(
            title='Save Log File',
            defaultextension='.log',
            filetypes=[
                ('Log files', '*.log'),
                ('Text files', '*.txt'),
                ('All files', '*.*'),
            ],
            initialfile=default_name,
            initialdir=self.work_dir,
        )
        if filepath:
            with open(filepath, 'w') as f:
                f.write(content)
                if not content.endswith('\n'):
                    f.write('\n')
            self.status_var.set(f'Log saved: {filepath}')

    # --- Window close ---

    def on_close(self):
        if self._process is not None:
            if not messagebox.askyesno('Program Running',
                    'A program is still running. Quit anyway?'):
                return
            self._process.terminate()
        if self.modified:
            if not messagebox.askyesno('Unsaved Changes',
                    'You have unsaved changes. Quit anyway?'):
                return
        self.root.destroy()


def main(templates_path=None):
    """Launch the Namelist Editor GUI."""
    if templates_path is None and len(sys.argv) > 1:
        templates_path = sys.argv[1]

    templates_dir = find_templates_dir(templates_path)

    if templates_dir is None:
        root = tk.Tk()
        root.withdraw()
        templates_dir = filedialog.askdirectory(
            title='Select EMsoftOO NamelistTemplates directory')
        root.destroy()
        if not templates_dir:
            print('Error: Could not find NamelistTemplates directory.')
            print('Usage: python -m emsoft.nml_editor /path/to/NamelistTemplates')
            sys.exit(1)

    root = tk.Tk()

    # High-DPI support on Windows
    if platform.system() == 'Windows':
        try:
            from ctypes import windll
            windll.shcore.SetProcessDpiAwareness(1)
        except Exception:
            pass

    app = NmlEditor(root, templates_dir)
    root.mainloop()


if __name__ == '__main__':
    main()
