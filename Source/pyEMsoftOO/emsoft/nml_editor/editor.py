"""
EMsoftOO Namelist Editor — a lightweight GUI for editing namelist template files.

Launch with:
    python -m emsoft.nml_editor
    python -m emsoft.nml_editor /path/to/NamelistTemplates
"""

import os
import sys
import platform
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


class NmlEditor:
    """Main application window."""

    def __init__(self, root, templates_dir):
        self.root = root
        self.templates_dir = templates_dir
        self.current_file = None
        self.modified = False
        self._highlight_job = None

        self.root.title('EMsoftOO Namelist Editor')
        self.root.geometry('820x650')
        self.root.minsize(600, 400)

        # Handle window close
        self.root.protocol('WM_DELETE_WINDOW', self.on_close)

        self._build_ui()
        self._bind_shortcuts()
        self._apply_theme()

    def _build_ui(self):
        """Create all UI elements."""
        # --- Top toolbar ---
        toolbar = ttk.Frame(self.root, padding=5)
        toolbar.pack(fill=tk.X)

        ttk.Label(toolbar, text='Program:').pack(side=tk.LEFT, padx=(0, 5))

        # Template selector
        self.templates = sorted(
            f for f in os.listdir(self.templates_dir) if f.endswith('.template')
        )
        display_names = [f.replace('.template', '') for f in self.templates]

        self.combo_var = tk.StringVar()
        self.combo = ttk.Combobox(toolbar, textvariable=self.combo_var,
                                  values=display_names, state='readonly', width=30)
        self.combo.pack(side=tk.LEFT, padx=(0, 10))
        self.combo.bind('<<ComboboxSelected>>', self.on_template_selected)

        ttk.Button(toolbar, text='Open .nml', command=self.open_file).pack(side=tk.LEFT, padx=2)
        ttk.Button(toolbar, text=f'Save .nml ({MOD_LABEL}+S)', command=self.save_file).pack(side=tk.LEFT, padx=2)

        # --- Editor area ---
        editor_frame = ttk.Frame(self.root)
        editor_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=(0, 0))

        # Line numbers
        self.linenums = tk.Text(editor_frame, width=4, padx=4, pady=4,
                                takefocus=0, border=0, state='disabled',
                                background='#f0f0f0', foreground='#999999')
        self.linenums.pack(side=tk.LEFT, fill=tk.Y)

        # Scrollbar
        scrollbar = ttk.Scrollbar(editor_frame)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Main text editor
        self.editor = tk.Text(editor_frame, wrap=tk.NONE, undo=True,
                              padx=6, pady=4, insertwidth=2,
                              yscrollcommand=self._on_scroll)
        self.editor.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self._scroll_both)

        # Configure editor font
        if IS_MAC:
            font = ('Menlo', 13)
        elif platform.system() == 'Windows':
            font = ('Consolas', 11)
        else:
            font = ('Monospace', 11)

        self.editor.configure(font=font)
        self.linenums.configure(font=font)

        # Track modifications
        self.editor.bind('<<Modified>>', self._on_modified)
        self.editor.bind('<KeyRelease>', self._schedule_highlight)

        # --- Status bar ---
        self.status_var = tk.StringVar(value='Select a program template to begin')
        status = ttk.Label(self.root, textvariable=self.status_var,
                           relief=tk.SUNKEN, anchor=tk.W, padding=(5, 2))
        status.pack(fill=tk.X, side=tk.BOTTOM)

        # --- Syntax highlighting tags ---
        self.editor.tag_configure('comment', foreground='#6a9955')
        self.editor.tag_configure('section', foreground='#2e8b57', font=font + ('bold',))
        self.editor.tag_configure('group', foreground='#0000cc', font=font + ('bold',))
        self.editor.tag_configure('param', foreground='#001080')
        self.editor.tag_configure('string', foreground='#a31515')
        self.editor.tag_configure('boolean', foreground='#7b2fa0')
        self.editor.tag_configure('number', foreground='#098658')

    def _bind_shortcuts(self):
        """Set up keyboard shortcuts."""
        mod = 'Command' if IS_MAC else 'Control'
        self.root.bind(f'<{mod}-s>', lambda e: self.save_file())
        self.root.bind(f'<{mod}-S>', lambda e: self.save_file())
        self.root.bind(f'<{mod}-o>', lambda e: self.open_file())
        self.root.bind(f'<{mod}-O>', lambda e: self.open_file())

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
        # Update scrollbar if it exists
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

        # Adjust width for large files
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

        # Remove existing tags
        for tag in ('comment', 'section', 'group', 'param', 'string',
                    'boolean', 'number'):
            self.editor.tag_remove(tag, '1.0', tk.END)

        content = self.editor.get('1.0', tk.END)
        lines = content.split('\n')

        for i, line in enumerate(lines):
            line_start = f'{i + 1}.0'
            line_end = f'{i + 1}.end'
            stripped = line.lstrip()

            if not stripped:
                continue

            # Namelist group start (&name)
            if stripped.startswith('&'):
                self.editor.tag_add('group', line_start, line_end)
                continue

            # Namelist terminator (/)
            if stripped == '/':
                self.editor.tag_add('group', line_start, line_end)
                continue

            # Section divider comments (!====)
            if stripped.startswith('!') and '===' in stripped:
                self.editor.tag_add('section', line_start, line_end)
                continue

            # Regular comments
            if stripped.startswith('!'):
                self.editor.tag_add('comment', line_start, line_end)
                continue

            # Lines with parameter = value
            if '=' in line and not stripped.startswith('!'):
                eq_pos = line.index('=')
                # Parameter name (before =)
                param_start = line_start
                param_end = f'{i + 1}.{eq_pos}'
                self.editor.tag_add('param', param_start, param_end)

                # Check for inline comment after the value
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

                # Highlight strings in value
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

                # Highlight booleans
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
            title += f' — {name}'
        if self.modified:
            title += ' (modified)'
        self.root.title(title)

    # --- Template loading ---

    def on_template_selected(self, event=None):
        if self.modified:
            if not messagebox.askyesno('Unsaved Changes',
                    'Current file has unsaved changes. Discard them?'):
                # Restore the combo to the previous selection
                return

        name = self.combo_var.get()
        template_file = name + '.template'
        filepath = os.path.join(self.templates_dir, template_file)

        if os.path.isfile(filepath):
            self._load_file(filepath)
            self.current_file = filepath
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

        # Remove trailing newline that tkinter adds
        if self.editor.get('end-2c', 'end-1c') == '\n':
            pass  # keep as-is

        self.editor.edit_modified(False)
        self.editor.edit_reset()  # clear undo stack
        self.modified = False
        self.current_file = filepath
        self._update_title()
        self._apply_highlighting()
        self.editor.see('1.0')

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
            initialdir=os.getcwd(),
        )
        if filepath:
            self._load_file(filepath)
            name = os.path.basename(filepath)
            line_count = int(self.editor.index('end-1c').split('.')[0])
            self.status_var.set(f'{name} loaded ({line_count} lines)')

    def save_file(self):
        """Save the editor content as a .nml file."""
        # Suggest a filename
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
            initialdir=os.getcwd(),
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
            self.status_var.set(f'Saved: {filepath}')

    # --- Window close ---

    def on_close(self):
        if self.modified:
            if not messagebox.askyesno('Unsaved Changes',
                    'You have unsaved changes. Quit anyway?'):
                return
        self.root.destroy()


def main(templates_path=None):
    """Launch the Namelist Editor GUI."""
    # Check command-line arguments
    if templates_path is None and len(sys.argv) > 1:
        templates_path = sys.argv[1]

    templates_dir = find_templates_dir(templates_path)

    if templates_dir is None:
        # Try to show a directory chooser
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
