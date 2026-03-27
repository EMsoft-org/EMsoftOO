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
        self.root.geometry('1000x700')
        self.root.minsize(700, 450)

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

        tk.Button(toolbar, text='  Open .nml  ', command=self.open_file).pack(side=tk.LEFT, padx=4)
        tk.Button(toolbar, text=f'  Save .nml ({MOD_LABEL}+S)  ',
                  command=self.save_file).pack(side=tk.LEFT, padx=4)

        # --- Main paned layout: template list on left, editor on right ---
        paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=(0, 0))

        # --- Left panel: template list with search ---
        left_frame = ttk.Frame(paned, padding=2)
        paned.add(left_frame, weight=0)

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

        # --- Right panel: editor ---
        right_frame = ttk.Frame(paned, padding=2)
        paned.add(right_frame, weight=1)

        editor_frame = ttk.Frame(right_frame)
        editor_frame.pack(fill=tk.BOTH, expand=True)

        # Line numbers
        self.linenums = tk.Text(editor_frame, width=4, padx=4, pady=4,
                                takefocus=0, border=0, state='disabled',
                                background='#f0f0f0', foreground='#999999',
                                font=self.editor_font)
        self.linenums.pack(side=tk.LEFT, fill=tk.Y)

        # Scrollbar
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

        # --- Status bar ---
        self.status_var = tk.StringVar(value='Select a template from the list to begin')
        status = ttk.Label(self.root, textvariable=self.status_var,
                           relief=tk.SUNKEN, anchor=tk.W, padding=(5, 2))
        status.pack(fill=tk.X, side=tk.BOTTOM)

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
                param_end = f'{i + 1}.{eq_pos}'
                self.editor.tag_add('param', line_start, param_end)

                # Find inline comment
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

                # Highlight strings
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
            title += f' \u2014 {name}'
        if self.modified:
            title += ' (modified)'
        self.root.title(title)

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
