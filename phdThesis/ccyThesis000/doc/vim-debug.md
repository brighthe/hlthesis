
The Terminal debugging plugin can be used to debug a program with gdb and view
the source code in a Vim window.  Since this is completely contained inside
Vim this also works remotely over an ssh connection.

When the |+terminal| feature is missing, the plugin will use the "prompt"
buffer type, if possible.  The running program will then use a newly opened
terminal window.  See |termdebug-prompt| below for details.


Starting ~
							*termdebug-starting*
Load the plugin with this command: >
	packadd termdebug
<							*:Termdebug*
To start debugging use `:Termdebug` or `:TermdebugCommand` followed by the
command name, for example: >
	:Termdebug vim

This opens two windows:

gdb window	A terminal window in which "gdb vim" is executed.  Here you
		can directly interact with gdb.  The buffer name is "!gdb".

program window	A terminal window for the executed program.  When "run" is
		used in gdb the program I/O will happen in this window, so
		that it does not interfere with controlling gdb.  The buffer
		name is "debugged program".

The current window is used to show the source code.  When gdb pauses the
source file location will be displayed, if possible.  A sign is used to
highlight the current position, using highlight group debugPC.

If the buffer in the current window is modified, another window will be opened
to display the current gdb position.  You can use `:Winbar` to add a window
toolbar there.

Focus the terminal of the executed program to interact with it.  This works
the same as any command running in a terminal window.

When the debugger ends, typically by typing "quit" in the gdb window, the two
opened windows are closed.

Only one debugger can be active at a time.
							*:TermdebugCommand*
If you want to give specific commands to the command being debugged, you can
use the `:TermdebugCommand` command followed by the command name and
additional parameters. >
	:TermdebugCommand vim --clean -c ':set nu'

Both the `:Termdebug` and `:TermdebugCommand` support an optional "!" bang
argument to start the command right away, without pausing at the gdb window
(and cursor will be in the debugged window).  For example: >
	:TermdebugCommand! vim --clean

To attach gdb to an already running executable or use a core file, pass extra
arguments.  E.g.: >
	:Termdebug vim core
	:Termdebug vim 98343

If no argument is given, you'll end up in a gdb window, in which you need to
specify which command to run using e.g. the gdb `file` command.


Example session ~
							*termdebug-example*
Start in the Vim "src" directory and build Vim: >
	% make
Make sure that debug symbols are present, usually that means that $CFLAGS
includes "-g".

Start Vim: >
	% ./vim

Load the termdebug plugin and start debugging Vim: >
	:packadd termdebug
	:Termdebug vim
You should now have three windows:
    source  - where you started, has a window toolbar with buttons
    gdb	    - you can type gdb commands here
    program - the executed program will use this window

You can use CTRL-W CTRL-W or the mouse to move focus between windows.
Put focus on the gdb window and type: >
	break ex_help
	run
Vim will start running in the program window. Put focus there and type: >
	:help gui
Gdb will run into the ex_help breakpoint.  The source window now shows the
ex_cmds.c file.  A red "1 " marker will appear in the signcolumn where the
breakpoint was set.  The line where the debugger stopped is highlighted.  You
can now step through the program.  Let's use the mouse: click on the "Next"
button in the window toolbar.  You will see the highlighting move as the
debugger executes a line of source code.

Click "Next" a few times until the for loop is highlighted.  Put the cursor on
the end of "eap->arg", then click "Eval" in the toolbar.  You will see this
displayed:
	"eap->arg": 0x555555e68855 "gui" ~
This way you can inspect the value of local variables.  You can also focus the
gdb window and use a "print" command, e.g.: >
	print *eap
If mouse pointer movements are working, Vim will also show a balloon when the
mouse rests on text that can be evaluated by gdb.

Now go back to the source window and put the cursor on the first line after
the for loop, then type: >
	:Break
You will see a ">>" marker appear, this indicates the new breakpoint.  Now
click "Cont" in the toolbar and the code until the breakpoint will be
executed.

You can type more advanced commands in the gdb window.  For example, type: >
	watch curbuf
Now click "Cont" in the toolbar (or type "cont" in the gdb window). Execution
will now continue until the value of "curbuf" changes, which is in do_ecmd().
To remove this watchpoint again type in the gdb window: >
	delete 3

You can see the stack by typing in the gdb window: >
	where
Move through the stack frames, e.g. with: >
	frame 3
The source window will show the code, at the point where the call was made to
a deeper level.


Stepping through code ~
							*termdebug-stepping*
Put focus on the gdb window to type commands there.  Some common ones are:
- CTRL-C	interrupt the program
- next		execute the current line and stop at the next line
- step		execute the current line and stop at the next statement,
		entering functions
- until		execute until past the current cursor line or past a specified
		position or the current stack frame returns
- finish	execute until leaving the current function
- where		show the stack
- frame N	go to the Nth stack frame
- continue	continue execution

						*:Run* *:Arguments*
In the window showing the source code these commands can be used to control
gdb:
 `:Run` [args]	    run the program with [args] or the previous arguments
 `:Arguments` {args}  set arguments for the next `:Run`

 *:Break*	set a breakpoint at the cursor position
 :Break {position}
		set a breakpoint at the specified position
 *:Clear*	delete the breakpoint at the cursor position

 *:Step*	execute the gdb "step" command
 *:Over*	execute the gdb "next" command (`:Next` is a Vim command)
 *:Until*	execute the gdb "until" command
 *:Finish*	execute the gdb "finish" command
 *:Continue*	execute the gdb "continue" command
 *:Stop*	interrupt the program

If 'mouse' is set the plugin adds a window toolbar with these entries:
  Step		`:Step`
  Next		`:Over`
  Finish	`:Finish`
  Cont		`:Continue`
  Stop		`:Stop`
  Eval		`:Evaluate`
This way you can use the mouse to perform the most common commands.  You need
to have the 'mouse' option set to enable mouse clicks.
								*:Winbar*
You can add the window toolbar in other windows you open with: >
  :Winbar

If gdb stops at a source line and there is no window currently showing the
source code, a new window will be created for the source code.  This also
happens if the buffer in the source code window has been modified and can't be
abandoned.

Gdb gives each breakpoint a number.  In Vim the number shows up in the sign
column, with a red background.  You can use these gdb commands:
- info break	list breakpoints
- delete N	delete breakpoint N
You can also use the `:Clear` command if the cursor is in the line with the
breakpoint, or use the "Clear breakpoint" right-click menu entry.


Inspecting variables ~
					*termdebug-variables* *:Evaluate*
 `:Evaluate`	    evaluate the expression under the cursor
 `K`		    same (see |termdebug_map_K| to disable)
 `:Evaluate` {expr}   evaluate {expr}
 `:'<,'>Evaluate`     evaluate the Visually selected text

This is similar to using "print" in the gdb window.
You can usually shorten `:Evaluate` to `:Ev`.


Other commands ~
							*termdebug-commands*
 *:Gdb*	     jump to the gdb window
 *:Program*    jump to the window with the running program
 *:Source*     jump to the window with the source code, create it if there
	     isn't one
 *:Asm*	     jump to the window with the disassembly, create it if there
	     isn't one

Events ~
							*termdebug-events*
Four autocommands can be used: >
	au User TermdebugStartPre  echomsg 'debugging starting'
	au User TermdebugStartPost echomsg 'debugging started'
	au User TermdebugStopPre   echomsg 'debugging stopping'
	au User TermdebugStopPost  echomsg 'debugging stopped'
<
						*TermdebugStartPre*
TermdebugStartPre		Before starting debugging.
				Not triggered if the debugger is already
				running or the debugger command cannot be
				executed.
						*TermdebugStartPost*
TermdebugStartPost		After debugging has initialized.
				If a "!" bang is passed to `:Termdebug` or
				`:TermdebugCommand` the event is triggered
				before running the provided command in gdb.
						*TermdebugStopPre*
TermdebugStopPre		Before debugging ends, when gdb is terminated,
				most likely after issuing a "quit" command in
				the gdb window.
						*TermdebugStopPost*
TermdebugStopPost		After debugging has ended, gdb-related windows
				are closed, debug buffers wiped out and
				the state before the debugging was restored.


Prompt mode ~
						*termdebug-prompt*
When the |+terminal| feature is not supported and on MS-Windows, gdb will run
in a buffer with 'buftype' set to "prompt".  This works slightly differently:
- The gdb window will be in Insert mode while typing commands.  Go to Normal
  mode with <Esc>, then you can move around in the buffer, copy/paste, etc.
  Go back to editing the gdb command with any command that starts Insert mode,
  such as `a` or `i`.
- The program being debugged will run in a separate window.  On MS-Windows
  this is a new console window.  On Unix, if the |+terminal| feature is
  available a Terminal window will be opened to run the debugged program in.

						*termdebug_use_prompt*
Prompt mode can be used even when the |+terminal| feature is present with: >
	let g:termdebug_config['use_prompt'] = 1
Or if there is no g:termdebug_config: >
	let g:termdebug_use_prompt = 1
<
						*termdebug_map_K*
The K key is normally mapped to :Evaluate. If you do not want this use: >
	let g:termdebug_config['map_K'] = 0
Or if there is no g:termdebug_config: >
	let g:termdebug_map_K = 0
<
						*termdebug_disasm_window*
If you want the Asm window shown by default, set the flag to 1.
the "disasm_window_height" entry can be used to set the window height: >
	let g:termdebug_config['disasm_window'] = 1
	let g:termdebug_config['disasm_window_height'] = 15
or, if there is no g:termdebug_config: >
	let g:termdebug_disasm_window = 15
Any value greater than 1 will set the Asm window height to that value.

Communication ~
						*termdebug-communication*
There is another, hidden, buffer, which is used for Vim to communicate with
gdb.  The buffer name is "gdb communication".  Do not delete this buffer, it
will break the debugger.

Gdb has some weird behavior, the plugin does its best to work around that.
For example, after typing "continue" in the gdb window a CTRL-C can be used to
interrupt the running program.  But after using the MI command
"-exec-continue"  pressing CTRL-C does not interrupt.  Therefore you will see
"continue" being used for the `:Continue` command, instead of using the
communication channel.


Customizing ~
				*termdebug-customizing* *g:termdebug_config*
In the past several global variables were used for configuration.  These are
deprecated, using the g:termdebug_config dictionary is preferred.  When
g:termdebug_config exists the other global variables will not be used.


GDB command ~
							*g:termdebugger*
To change the name of the gdb command, set "debugger" entry in
g:termdebug_config or the "g:termdebugger" variable before invoking
`:Termdebug`: >
	let g:termdebug_config['command'] = "mygdb"
Or if there is no g:termdebug_config: >
	let g:termdebugger = "mygdb"

If the command needs an argument use a List: >
	let g:termdebug_config['command'] = ['rr', 'replay', '--']
Or if there is no g:termdebug_config: >
	let g:termdebugger = ['rr', 'replay', '--']

Several arguments will be added to make gdb work well for the debugger.
If you want to modify them, add a function to filter the argument list: >
	let g:termdebug_config['command_filter'] = MyDebugFilter

If you do not want the arguments to be added, but you do need to set the
"pty", use a function to add the necessary arguments: >
	let g:termdebug_config['command_add_args'] = MyAddArguments
The function will be called with the list of arguments so far, and a second
argument that is the name of the pty.
							*gdb-version*
Only debuggers fully compatible with gdb will work.  Vim uses the GDB/MI
interface.  The "new-ui" command  requires gdb version 7.12 or later.  if you
get this error:
	Undefined command: "new-ui". Try "help".~
Then your gdb is too old.


Colors ~
						*hl-debugPC* *hl-debugBreakpoint*
The color of the signs can be adjusted with these highlight groups:
- debugPC		the current position
- debugBreakpoint	a breakpoint

The defaults are, when 'background' is "light":
  hi debugPC term=reverse ctermbg=lightblue guibg=lightblue
  hi debugBreakpoint term=reverse ctermbg=red guibg=red

When 'background' is "dark":
  hi debugPC term=reverse ctermbg=darkblue guibg=darkblue
  hi debugBreakpoint term=reverse ctermbg=red guibg=red


Shortcuts ~
							*termdebug_shortcuts*

You can define your own shortcuts (mappings) to control gdb, that can work in
any window, using the TermDebugSendCommand() function.  Example: >
	map ,w :call TermDebugSendCommand('where')<CR>
The argument is the gdb command.


Popup menu ~
							*termdebug_popup*

By default the Termdebug plugin sets 'mousemodel' to "popup_setpos" and adds
these entries to the popup menu:
	Set breakpoint		`:Break`
	Clear breakpoint	`:Clear`
	Evaluate		`:Evaluate`
If you don't want this then disable it with: >
	let g:termdebug_config['popup'] = 0
or if there is no g:termdebug_config: >
	let g:termdebug_popup = 0


Vim window width ~
							*termdebug_wide*

To change the width of the Vim window when debugging starts and use a vertical
split: >
	let g:termdebug_config['wide'] = 163
Or if there is no g:termdebug_config: >
	let g:termdebug_wide = 163

This will set 'columns' to 163 when `:Termdebug` is used.  The value is
restored when quitting the debugger.

If the wide value is set and 'columns' is already a greater value, then a
vertical split will be used without modifying 'columns'.

Set the wide value to 1 to use a vertical split without ever changing
'columns'.  This is useful when the terminal can't be resized by Vim.


 vim:tw=78:ts=8:noet:ft=help:norl:

