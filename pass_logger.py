import getopt
import sys
import tkinter as tk
from tkinter import *
from tkinter import filedialog
from tkinter.scrolledtext import ScrolledText
from unicodedata import lookup
import os

class Diacritical:
    """Mix-in class that adds keyboard bindings for accented characters, plus
    other common functionality.

    An inheriting class must define a select_all method that will respond
    to Ctrl-A."""

    accents = (('acute', "'"), ('grave', '`'), ('circumflex', '^'),
               ('tilde', '='), ('diaeresis', '"'), ('cedilla', ','),
               ('stroke', '/'), ('ring above', ';'))

    def __init__(self):
        # Fix some key bindings
        self.bind("<Control-Key-a>", self.select_all)
        
        # We will need Ctrl-/ for the "stroke", but it cannot be unbound, so
        # let's prevent it from being passed to the standard handler
        #self.bind("<Control-Key-/>", lambda event: "break")
        # Diacritical bindings
        for a, k in self.accents:
            # Little-known feature of Tk, it allows to bind an event to
            # multiple keystrokes
            #self.bind("<Control-Key-%s><Key>" % k,
            #          lambda event, a=a: self.insert_accented(event.char, a))
            0

    def insert_accented(self, c, accent):
        if c.isalpha():
            if c.isupper():
                cap = 'capital'
            else:
                cap = 'small'
            try:
                c = lookup("latin %s letter %c with %s" % (cap, c, accent))
                self.insert(INSERT, c)
                # Prevent plain letter from being inserted too, tell Tk to
                # stop handling this event
                return "break"
            except KeyError as e:
                pass

class DiacriticalEntry(Entry, Diacritical):
    """Tkinter Entry widget with some extra key bindings for
    entering typical Unicode characters - with umlauts, accents, etc."""

    def __init__(self, master=None, **kwargs):
        Entry.__init__(self, master, **kwargs)
        Diacritical.__init__(self)

    def select_all(self, event=None):
        self.selection_range(0, END)
        return "break"

class DiacriticalText(ScrolledText, Diacritical):
    """Tkinter ScrolledText widget with some extra key bindings for
    entering typical Unicode characters - with umlauts, accents, etc."""

    def __init__(self, master=None, **kwargs):
        ScrolledText.__init__(self, master, **kwargs)
        Diacritical.__init__(self)

    def select_all(self, event=None):
        self.tag_add(SEL, "1.0", "end-1c")
        self.mark_set(INSERT, "1.0")
        self.see(INSERT)
        return "break"

class CombinedWidget:
    def __init__(self, master, inc_list, key, counter):
        self.inc_list = inc_list
        self.index_counter = 0
        self.key = key
        self.label = tk.Label(master, text=key + "[0/%d]" % (len(inc_list)-1))
        self.label.grid()
        self.entry = DiacriticalText(master, width=50, height=10)
        try:
            self.entry.insert(tk.END, "\n".join(inc_list[self.index_counter]))
        except IndexError:
            self.entry.insert(tk.END, "@@@@@@No results found!@@@@@@")
        self.entry.grid()
        self.button = tk.Button(master, text="Next", command=self.next_index)
        self.button.grid()
        self.button_2 = tk.Button(master, text="Previous", command=self.previous_index)
        self.button_2.grid()

    def next_index(self, event=None):
        self.index_counter += 1
        if self.index_counter >= len(self.inc_list):
            self.index_counter = 0
        self.entry.delete('1.0', tk.END)
        try:
            self.entry.insert(tk.END, "\n".join(self.inc_list[self.index_counter]))
        except IndexError:
            self.entry.insert(tk.END, "@@@@@@No results found!@@@@@@")
        self.label.config(text = self.key + "[%d/%d]" % (self.index_counter, len(self.inc_list)-1))
    def previous_index(self, event=None):
        self.index_counter -= 1
        if self.index_counter < 0:
            self.index_counter = len(self.inc_list) - 1
        self.entry.delete('1.0', tk.END)
        try:
            self.entry.insert(tk.END, "\n".join(self.inc_list[self.index_counter]))
        except IndexError:
            self.entry.insert(tk.END, "@@@@@@No results found!@@@@@@")
        self.label.config(text = self.key + "[%d/%d]" % (self.index_counter, len(self.inc_list)-1))

def find_first_last(l):
    first = ""
    last = ""
    mark = False

    for i in range(len(l)):
        if mark:
            break
        error = False
        for k in range(len(l[i])):
            if "returned error" in l[i][k]:
                error = True
                break
        if error:
            error = False
            continue
        else:
            mark = True
            first = l[i]
            break
    
    mark = False

    for i in range(len(l)):
        if mark:
            break
        error = False
        reversed_index = len(l)-(1+i)
        for k in range(len(l[reversed_index])):
            if "returned error" in l[reversed_index][k]:
                error = True
                break
        if error:
            error = False
            continue
        else:
            mark = True
            last = l[reversed_index]
    
    return first, last

def combined_string(result_dict, key):
    if key == 'cubeadcs combined':
        combined_cmd_list = [
            'cubeadcs tlm',
            'cubeadcs get gyro_rates',
            'cubeadcs get attitude',
            'cubeadcs get estimated_rates',
            'cubeadcs cw get',
            'cubeadcs get acp_state'
            ]
    elif key == 'ttnc_hk combined':
        combined_cmd_list = [
                'ax100 hk',
                'ax2150 hk'
            ]
    elif key == 'obc_hk combined':
        combined_cmd_list = [
                'rparam init 1 4',
                'rparam init 1 1'
            ]
    else:
        raise Exception("Unknown key to combined! %s" % key)
    
    result_list = []
    for j in range(len(result_dict[combined_cmd_list[0]])):
        combined_string = []
        for i in range(len(combined_cmd_list)):
            found = False
            current_list = result_dict[combined_cmd_list[i]]
            #print(current_list)
            try:
                if "returned error" in "\n".join(current_list[j]):
                    for k in range(len(current_list)):
                        if "returned error" in "\n".join(current_list[k]):
                            continue
                        else:
                            found = True
                            combined_string.append("\n".join(current_list[k]))
                    if found is False:
                        combined_string.append("\n".join(current_list[k]))
                else:
                    #raise Exception
                    combined_string.append("\n".join(current_list[j]))
            except IndexError as e:
                #print(e)
                continue
        result_list.append(combined_string)
    
    # remove duplicates

    for i in range(len(result_list)):
        if i == len(result_list) - 1:
            break
        else:
            try:
                if "\n".join(result_list[i]) == "\n".join(result_list[i+1]):
                    result_list.pop(i)
            except IndexError:
                continue

    return result_list




def parse_dialog(first_beacon, last_eps_hk, cubeadcs_combined, obc_hk, uptime):
    # find slgv2 current draw from first beacon. 
    slgv2_aos_current_draw = None
    for i in first_beacon:
        if "(H1-47)" in i:
            slgv2_aos_current_draw = int(i.split("[")[-1].split(",")[0].strip())

    #print(last_eps_hk)
    battery_voltage = None
    slgv2_LOS_current_draw = None
    avg_temp = None
    for i in last_eps_hk:
        if "mV" in i and "EN" not in i:
            battery_voltage = i.split("|")[1].strip()
        if "(H1-47)" in i:
            slgv2_LOS_current_draw = int(i.split("[")[-1].split(",")[0].strip())
        if "Temp" in i:
            avg_temp = i.split(" ")[-1]
    
    #print(cubeadcs_combined[0].split("\n"))
    tlm = cubeadcs_combined[0].split("\n")
    acp_run_mode = None
    estimation_mode = None
    control_mode = None
    com_temp = None
    for i in tlm:
        if "ACP run mode" in i:
            #print(i.split(":")[-1].strip())
            acp_run_mode = i.split(":")[-1].strip()
        if "Current Estimation mode" in i:
            #print(i.split(":")[-1].strip())
            estimation_mode = i.split(":")[-1].strip()
        if "Current Control mode" in i:
            #print(i.split(":")[-1].strip())
            control_mode = i.split(":")[-1].strip()
        if "CubeComputer temperature" in i:
            #print(i.split("=")[-1].strip())
            com_temp = i.split("=")[-1].strip()
    
    #print(cubeadcs_combined[1].split("\n"))
    gyro = cubeadcs_combined[1].split("\n")
    xrate = None
    yrate = None
    zrate = None
    for i in gyro:
        if "Xrate" in i:
            xrate = i.split(":")[-1].strip()
        if "Yrate" in i:
            yrate = i.split(":")[-1].strip()
        if "Zrate" in i:
            zrate = i.split(":")[-1].strip()
    
    #print(cubeadcs_combined[2].split("\n"))
    attitude = cubeadcs_combined[2].split("\n")
    est_roll = None
    est_pitch = None
    est_yaw = None
    for i in attitude:
        if "Estimated Roll" in i:
            est_roll = i.split(":")[-1].strip()
        if "Estimated Pitch" in i:
            est_pitch = i.split(":")[-1].strip()
        if "Estimated Yaw" in i:
            est_yaw = i.split(":")[-1].strip()
    

    #print(cubeadcs_combined[3].split("\n"))
    est_gyro = cubeadcs_combined[3].split("\n")
    est_xrate = None
    est_yrate = None
    est_zrate = None
    for i in est_gyro:
        if "Xrate" in i:
            est_xrate = i.split(":")[-1].strip()
        if "Yrate" in i:
            est_yrate = i.split(":")[-1].strip()
        if "Zrate" in i:
            est_zrate = i.split(":")[-1].strip()
    

    #print(obc_hk[0].split("\n"))
    rparam_1_4 = obc_hk[0].split("\n")
    ram_image = None
    bootcount = None
    for i in rparam_1_4:
        if "ram_image" in i:
            ram_image = i.split()[-1]
        if 'bootcount' in i:

            bootcount = i.split()[-1]

    #print(obc_hk[1].split("\n"))
    rparam_1_1 = obc_hk[1].split("\n")
    swload_cnt1 = None
    for i in rparam_1_1:
        if "swload_cnt1" in i:
            swload_cnt1 = i.split()[-1]
    
    uptime_value = uptime[-1].split("is")[-1].strip()

    line_1 = "SLGv2 current draw at AOS: %smA (first beacon)\n" % slgv2_aos_current_draw
    line_2 = "SLGv2 last known current draw prior to LOS: %smA (last eps hk)\n" % slgv2_LOS_current_draw
    line_3 = "Battery Normal (%s), %s average EPS temperature\n" % (battery_voltage, avg_temp)
    line_4 = "tumbling rates============\nGyro: %s, %s, %s\n" % (xrate, yrate, zrate)
    line_5 = "Estimated: %s, %s, %s\n" % (est_xrate, est_yrate, est_zrate)
    line_6 = "Roll: %s, Pitch: %s, Yaw: %s\n" % (est_roll, est_pitch, est_yaw)
    line_7 = "ADCS modes: run = %s, estimation = %s, control = %s, CubeComptuerTemp = %s\n" % (acp_run_mode, estimation_mode, control_mode, com_temp)
    line_8 = "OBC ram_image = %s, OBC bootcount = %s, swload_cnt1=%s, uptime = %s\n" % (ram_image, str(bootcount), str(swload_cnt1), uptime_value)

    combined_line= [line_1 , line_2 , line_3 , line_4 , line_5 , line_6 , line_7, line_8]
    #print(combined_line)
    return [combined_line]

def main(argv):
    window = tk.Tk()

    def donothing():
        x = 0

    def parse_file():
        filename = filedialog.askopenfilename(initialdir = "%s" % os.getcwd(),
                                          title = "Select a File",
                                          filetypes = (("Text files",
                                                        "*.txt*"),
                                                       ("all files",
                                                        "*.*")))
      
        # Change label contents
        with open(filename, 'r') as f:
            txt = f.read()

        commands = []

        result_dict = {
            "Beacon": [],
            "eps hk": [],
            "ax100 hk": [],
            "ax2150 hk": [],
            "rparam init 1 4": [],
            "rparam init 1 1": [],
            "uptime 1": [],
            "obc timesync": [],
            "cubeadcs tlm": [],
            "cubeadcs get gyro_rates": [],
            "cubeadcs get attitude": [],
            "cubeadcs get estimated_rates": [],
            "cubeadcs cw get": [],
            "cubeadcs get acp_state": [],
            "slgv2 status print": [],
            "cubeadcs combined": [],
            "ttnc_hk combined": [],
            "obc_hk combined": [],
            }

        txt = txt.split("\n")
        for i in range(len(txt)):
            line = txt[i]
            intermediate = []
            if "csp-term #" in line and len(line) > 11:
                # means that a legit command was found.
                try: 
                    if "csp: OUT" in txt[i+1]:
                        command_unix = float(txt[i+1].split("I")[0])
                except:
                    0
                intermediate.append(line)
                # we want to scan the results of the command until the next csp-term #
                error = False
                for k in range(1, 100):
                    try:
                        if "returned error" in txt[i+k]:
                            error = True
                        if "csp-term #" in txt[i+k]:
                            if "rparam getall" in txt[i+k]:
                                continue
                            if "uptime 1" in txt[i+k]:
                                continue
                            if "obc timesync" in txt[i+k]:
                                continue
                            else:
                                i = i+k-1
                                break
                        else:
                            intermediate.append(txt[i+k])
                    except:
                        0
                
                add_on = ""
                if "ZMQHUB" not in line:
                    if error:
                        add_on = "(RS ERROR/TIMEOUT)"
                    commands.append(line.replace("csp-term #", "") + add_on)

                # now we got the intermediate, we go through and switch case based on the commands.

                for key in result_dict.keys():
                    if key in line:
                        result_dict[key].append(intermediate)
                    else:
                        # unknown command?
                        # print(line)
                        0

            elif "Beacon" in line:
                for k in range(1, 100):
                    if "bcn-handler: csv" in txt[i+k]:
                        intermediate.append(txt[i+k])
                        i = i+k-1
                        break
                    else:
                        intermediate.append(txt[i+k])
                
                result_dict["Beacon"].append(intermediate)
            
            else:
                for key in result_dict.keys():
                    if key in line:
                        
                        intermediate.append(line)
                        # we want to scan the results of the command until the next csp-term #
                        for k in range(1, 100):
                            if "csp-term #" in txt[i+k]:
                                if "rparam getall" in txt[i+k]:
                                    continue
                                else:
                                    i = i+k-1
                                    break
                            else:
                                intermediate.append(txt[i+k])
                        result_dict[key].append(intermediate)
                        break

        
        result_dict['cubeadcs combined'] = combined_string(result_dict, 'cubeadcs combined')
        result_dict['ttnc_hk combined'] = combined_string(result_dict, 'ttnc_hk combined')
        result_dict['obc_hk combined'] = combined_string(result_dict, 'obc_hk combined')

        first_eps_hk, last_eps_hk = find_first_last(result_dict['eps hk'])
    
        try:
            result_dict['parse_dialog'] = parse_dialog(result_dict['Beacon'][0], last_eps_hk, result_dict['cubeadcs combined'][0], result_dict['obc_hk combined'][0], result_dict['uptime 1'][0])
        except:
            0
        
                
        counter = 0
        bg_list = ['cyan', 'gray2', 'white', 'lavender']
        for key in result_dict.keys():
            if "combined" in key:
                background = 'cyan'
            else:
                background = 'lavender'
            frame = tk.Frame(window, bg=background, pady=4)
            frame.grid(row = counter%4, column=int(counter/4))
            new = CombinedWidget(frame, result_dict[key], key, counter)
            counter += 1
        

        frame = tk.Frame(window, bg=bg_list[counter%4], pady=4)
        frame.grid(row = counter%4, column=int(counter/4))
        label = tk.Label(frame, text="Commands Sent")
        label.grid()
        entry = DiacriticalText(frame, width=30, height=20)
        entry.insert(tk.END, "\n".join(commands))
        entry.grid()
        counter += 1


    def help_win():
        help = tk.Toplevel()
        help.title("How to Use")
        about= """
        1. Click File -> Open and select the txt. Should be in NuLIoN Telems
        2. Multiple windows will be loaded. 
        3. Title of each window is <COMMAND> and [CURRENT_INDEX/TOTAL_INDEX]
        4. A value of TOTAL_INDEX = -1 means nothing found. 
        4. Click Next and Previous to cycle through the index.
        5. Click into the widget
        6. and use control+A to select all and then control+C to copy.
        """
        t = DiacriticalText(help)
        t.insert(tk.END, about)
        t.pack()
        tk.Button(help, text="OK", command=help.destroy).pack()
    
    def about_win():
        about = tk.Toplevel()
        about.title("FRAMEWORK")
        entry = """
        _   _           ____                                
        | \ | |  _   _  / ___|   _ __     __ _    ___    ___ 
        |  \| | | | | | \___ \  | '_ \   / _` |  / __|  / _ \
        | |\  | | |_| |  ___) | | |_) | | (_| | | (__  |  __/
        |_| \_|  \__,_| |____/  | .__/   \__,_|  \___|  \___|
                                |_|                          
        (C)2024 NuSpace Pte Ltd

        Description:
            Parse that beacon!

        Author: Hubert Khoo Hui Bo
        Date Written: 3rd May 2024
        """
        t = DiacriticalText(about, width=200)
        t.insert(tk.END, entry)
        t.pack()
        tk.Button(about, text="OK", command=about.destroy).pack()

    menubar = Menu(window)
    filemenu = Menu(menubar, tearoff=0)
    filemenu.add_command(label="Open", command=parse_file)
    filemenu.add_separator()
    filemenu.add_command(label="Exit", command=window.quit)
    menubar.add_cascade(label="File", menu=filemenu)

    helpmenu = Menu(menubar, tearoff=0)
    helpmenu.add_command(label="How to Use", command=help_win)
    helpmenu.add_command(label="About...", command=about_win)
    menubar.add_cascade(label="Help", menu=helpmenu)

    window.config(menu=menubar)

    window.mainloop()
    
if __name__ == "__main__":
    main(sys.argv[1:])