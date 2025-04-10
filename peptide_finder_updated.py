from Bio import SeqIO
import ast
import os
import sys
from datetime import date,timedelta
from random import randint
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import tkinter.scrolledtext as st 
import threading
from idlelib.tooltip import Hovertip
import time  
today = str(date.today())

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("TPM")
        self.iconbitmap(resource_path("helix-removebg-preview.ico"))
        self.geometry("255x255")
        self.resizable(False,False)
        self.selected_file = ''
        self.create_widgets()

    def destroy_all(self):
        for widget in self.winfo_children():
            widget.destroy()
    def update_selected_file(self):
        starting = self.selected_file[self.selected_file.index('sequences/')+10:]
        display_file = ''
        per_line = 34
        for i in range(0,len(starting),per_line):
            if ((i+per_line) > len(starting)):
                display_file += starting[i:]
            else:
                display_file += starting[i:i+per_line]+'\n'
        if (len(display_file) <=33):
            display_file += '\n'
        if len(display_file)>=65:
            display_file = display_file[0:64]+'[...]'
        return display_file
            
    def error_reset(self):
        self.destroy_all()
        self.title('TPM: Error')
        self.geometry('255x255')
        self.error_label = tk.Label(self, text="Error encountered!\nPlease make sure you are\nusing valid parameters.",
                                    font=("Helvetica",11),fg='red')
        self.error_label.place(relx=0.5,rely=0.3,anchor='center')
        self.error_home_btn = tk.Button(self, text="Return Home", command=self.create_widgets)
        self.error_home_btn.place(relx=0.5,rely=0.5,anchor='center')
        
    def browseFiles(self):
        self.selected_file = tk.filedialog.askopenfilename(initialdir = os.getcwd()+'/sequences',
                                              title = "Select a File",
                                              filetypes = (("FASTA files",
                                                            "*.fasta*"),("Text files",
                                                            "*.txt*")
                                                           ))
        self.create_widgets()

    def validateDirectories(self):
        return (os.path.isdir('sequences') and os.path.isdir('outputs') and
                os.path.isdir('user_data') and os.path.exists('user_data/protein_profiles.txt')
                and os.path.exists('user_data/peptide_profiles.txt'))
    def createDir(self):
        for path in ('sequences','outputs','user_data'):
            if not os.path.isdir(path):
                os.mkdir(path)
        for file in ('user_data/peptide_profiles.txt','user_data/protein_profiles.txt'):
            if not os.path.exists(file):
                with open(file,'w') as z:
                    pass
        self.create_widgets()
    def buildDirectories(self):
        self.destroy_all()
        self.title('TPM: Setup')
        self.top_label = tk.Label(self, text="Setup Required",font=("Helvetica",10,"bold"))
        self.top_label.pack(pady=2)
        self.mid_label = tk.Label(self, text="TPM needs to create local directories\nand files to store user parameters,\nsequence files, and output files.")
        self.mid_label.pack(pady=2)
        self.build_btn = tk.Button(self, text="Finish Setup", command=self.createDir)
        self.build_btn.pack(pady=2)
        
    
                
    def create_widgets(self):
        self.destroy_all()
        self.title("TPM")
        self.geometry("255x255")
        self.top_label = tk.Label(self, text="Select sequences for analysis:",font=("Helvetica",10,"bold"))
        self.top_label.pack(pady=2)
        self.start_button = tk.Button(self, text="Browse Files", command=self.browseFiles)
        self.start_button.pack(pady=2)
        self.sel_lbl = tk.Label(self, text="Selected File:",font = ("Helvetica", 9, "bold"))
        self.sel_lbl.pack(pady=(2,0))
        self.selected_file_lbl = tk.Label(self, text="[none]\n",anchor='n')
        if not (len(self.selected_file)==0):
            self.selected_file_lbl.configure(text=self.update_selected_file(),font = ("Courier", 9, "bold"),anchor='n')
        self.selected_file_lbl.pack(pady=(0,3))
        self.genome_to_prot_btn = tk.Button(self, text="Genome to Target Protein", command=self.proteinTranslate)
        self.genome_to_prot_btn.pack(pady=2)
        self.prot_to_pep_btn = tk.Button(self, text="Protein to Peptide Report", command=self.peptideReport)
        self.prot_to_pep_btn.pack(pady=2)
        if (len(self.selected_file)==0):
            self.prot_to_pep_btn['state'] = 'disabled'
            self.genome_to_prot_btn['state'] = 'disabled'

        ve2_frame = tk.Frame(self,width=255,height=16)
        ve2_frame.pack(pady=1,anchor='s')
        self.viewEdit_label = tk.Label(self, text="View/Edit User Sequences",font=("Helvetica",9,"bold"))
        self.viewEdit_label.pack(pady=(4,1))
        ve_frame = tk.Frame(self,width=255,height=25)
        ve_frame.pack(pady=1,anchor='s')
        self.edit_prot_btn = tk.Button(ve_frame, text="Protein", command=self.editProteinsProfiles)
        self.edit_prot_btn.grid(row=0,column=0,padx=3)
        self.edit_pep_btn = tk.Button(ve_frame, text="Peptide", command=self.editPeptidesProfiles)
        self.edit_pep_btn.grid(row=0,column=1,padx=3)
        
        if not self.validateDirectories():
            self.buildDirectories()

    def editProteinsProfiles(self):
        self.protein_profile_dict = {}
        self.protein_profile_dict['None'] = ''
        try:
            with open('user_data/protein_profiles.txt') as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    self.protein_profile_dict[record.id] = record.seq
        except:
            pass
        self.protKeys = list(self.protein_profile_dict.keys())
        self.title("TPM: View/Edit")
        self.geometry('255x165')
        self.destroy_all()
        self.title_lbl = tk.Label(self, text="Saved Proteins",font=("Helvetica",10,"bold"))
        self.title_lbl.grid(row=0,column=0,columnspan=2)
        self.return_home_btn = tk.Button(self, text="\u2BAA"+" Back", command=self.create_widgets)
        self.return_home_btn.grid(row=0,column=0,sticky='w')
        self.columnconfigure(0,weight=1)
        self.columnconfigure(1,weight=4)  
        self.prot_label = tk.Label(self, text="Protein\nSequence:")
        self.prot_label.grid(row=1,column=0)
        self.prot_entry = st.ScrolledText(self,width=23,height=5,font=("Courier",9))
        self.prot_entry.grid(row=1,column=1)
        prot_tip = Hovertip(self.prot_label,'Select sequence to view/edit from the,\ndropdown menu, or create a new sequence.') 
        self.chosen_prot = tk.StringVar(self)
        self.prot_options = ttk.OptionMenu(self,self.chosen_prot, 'None', *self.protKeys, command = self.ve_prot_change)
        self.prot_options.grid(row=2,column=1)
        prot_bottom_frame = tk.Frame(self,width=255,height=25)
        prot_bottom_frame.grid(row=3,column=0,columnspan=2,pady=5,sticky='s')
        self.del_prot_btn = tk.Button(prot_bottom_frame, text="Delete", command=self.del_cur_prot)
        self.del_prot_btn.grid(row=0,column=0,padx=3)
        del_prot_tip = Hovertip(self.del_prot_btn,'Delete currently selected\nprotein sequence.')
        self.save_prot_btn = tk.Button(prot_bottom_frame, text="Save", command=self.edit_cur_prot)
        self.save_prot_btn.grid(row=0,column=1,padx=3)
        save_prot_tip = Hovertip(self.save_prot_btn,'Save changes to currently\nselected protein sequence.')
        self.new_prot_btn = tk.Button(prot_bottom_frame, text="New", command=self.new_prot_prof)
        self.new_prot_btn.grid(row=0,column=3,padx=3)
        new_prot_tip = Hovertip(self.new_prot_btn,'Write a new protein sequence.')
    def del_cur_prot(self):
        if self.chosen_prot.get() == 'None':
            return
        full_profs = ''
        with open('user_data/protein_profiles.txt','w') as f:
            for prot in self.protein_profile_dict:
                if prot != 'None' and prot != self.chosen_prot.get():
                    full_profs += '>'+prot+'\n'+self.protein_profile_dict[prot]+'\n'
            print(full_profs,file=f)
            del self.protein_profile_dict[self.chosen_prot.get()]
            self.protKeys = list(self.protein_profile_dict.keys())
            self.prot_entry.delete('1.0',tk.END)
            self.prot_options.destroy()
            self.prot_options = ttk.OptionMenu(self,self.chosen_prot, 'None', *self.protKeys, command = self.ve_prot_change)
            self.prot_options.grid(row=2,column=1)
    def edit_cur_prot(self):
        if self.chosen_prot.get() == 'None':
            return
        full_profs = ''
        with open('user_data/protein_profiles.txt','w') as f:
            for prot in self.protein_profile_dict:
                if prot == self.chosen_prot.get():
                    full_profs += '>'+prot+'\n'+self.prot_entry.get('1.0',tk.END).upper().replace('\n','').replace(' ','')+'\n'
                elif prot != 'None':
                    full_profs += '>'+prot+'\n'+self.protein_profile_dict[prot]+'\n'
            print(full_profs,file=f)
            self.protein_profile_dict[self.chosen_prot.get()] = self.prot_entry.get('1.0',tk.END)
            self.protKeys = list(self.protein_profile_dict.keys())
            self.prot_options.destroy()
            self.prot_options = ttk.OptionMenu(self,self.chosen_prot, self.chosen_prot.get(), *self.protKeys, command = self.ve_prot_change)
            self.prot_options.grid(row=2,column=1)
    def new_prot_prof(self):
        self.destroy_all()
        self.title_lbl = tk.Label(self, text="New Protein",font=("Helvetica",10,"bold"))
        self.title_lbl.grid(row=0,column=0,columnspan=2)
        self.return_home_btn = tk.Button(self, text="\u2BAA"+" Back", command=self.editProteinsProfiles)
        self.return_home_btn.grid(row=0,column=0,sticky='w')
        self.columnconfigure(0,weight=1)
        self.columnconfigure(1,weight=4)
        self.npname_lbl = tk.Label(self, text="Protein Name:")
        self.npname_lbl.grid(row=1,column=0)
        self.npname_entry = tk.Entry(self,width=23)
        self.npname_entry.grid(row=1,column=1)
        self.newprot_label = tk.Label(self, text="Protein\nSequence:")
        self.newprot_label.grid(row=2,column=0)
        newprotname_tip = Hovertip(self.npname_lbl,'Enter the name for the new protein.')
        self.newprot_entry = st.ScrolledText(self,width=23,height=5,font=("Courier",9))
        self.newprot_entry.grid(row=2,column=1)
        newprot_tip = Hovertip(self.newprot_label,'Type or paste in a new protein sequence.')
        self.save_new_prot_btn = tk.Button(self, text="Save", command=self.saveNewProtProf)
        self.save_new_prot_btn.grid(row=3,column=1)
    def saveNewProtProf(self):
        self.protein_profile_dict[self.npname_entry.get().replace(' ','_')] = self.newprot_entry.get('1.0',tk.END).upper().replace('\n','').replace(' ','')
        self.protKeys = list(self.protein_profile_dict.keys())
        with open('user_data/protein_profiles.txt','w') as f:
            for prot in self.protein_profile_dict:
                print('>'+prot+'\n'+self.protein_profile_dict[prot]+'\n',file=f)
        self.save_new_prot_btn.configure(text='Saved',state=tk.DISABLED)
    def ve_prot_change(self,*args):
        self.prot_entry.delete('1.0',tk.END)
        self.prot_entry.insert(tk.END, self.protein_profile_dict[self.chosen_prot.get()])

    def editPeptidesProfiles(self):
        self.title("TPM: View/Edit")
        self.destroy_all()
        self.geometry('255x165')
        self.title_lbl = tk.Label(self, text="Saved Peptides",font=("Helvetica",10,"bold"))
        self.title_lbl.grid(row=0,column=0,columnspan=2)
        self.return_home_btn = tk.Button(self, text="\u2BAA"+" Back", command=self.create_widgets)
        self.return_home_btn.grid(row=0,column=0,sticky='w')
        self.columnconfigure(0,weight=1)
        self.columnconfigure(1,weight=4)
        self.peptide_profile_dict = {}
        try:
            with open('user_data/peptide_profiles.txt','r') as file:
                self.peptide_profile_dict = ast.literal_eval(file.read())
        except:
            pass
        self.peptide_profile_dict['None'] = ''
        self.pepKeys = list(self.peptide_profile_dict.keys())
        self.pep_label = tk.Label(self, text="Peptide\nSequences:")
        self.pep_label.grid(row=1,column=0)
        self.pep_entry = st.ScrolledText(self,width=23,height=5,font=("Courier",9))
        self.pep_entry.grid(row=1,column=1)
        pep_tip = Hovertip(self.pep_label,'Select a peptide list to view/edit from the\ndropdown menu, or create a new one.') 
        self.chosen_pep = tk.StringVar(self)
        self.pep_options = ttk.OptionMenu(self,self.chosen_pep, 'None', *self.pepKeys, command = self.ve_pep_change)
        self.pep_options.grid(row=2,column=1)
        pep_bottom_frame = tk.Frame(self,width=255,height=25)
        pep_bottom_frame.grid(row=3,column=0,columnspan=2,pady=5,sticky='s')
        self.del_pep_btn = tk.Button(pep_bottom_frame, text="Delete", command=self.del_cur_pep)
        self.del_pep_btn.grid(row=0,column=0,padx=3)
        del_pep_tip = Hovertip(self.del_pep_btn,'Delete currently selected\npeptide list.')
        self.save_pep_btn = tk.Button(pep_bottom_frame, text="Save", command=self.edit_cur_pep)
        self.save_pep_btn.grid(row=0,column=1,padx=3)
        save_pep_tip = Hovertip(self.save_pep_btn,'Save changes to currently\nselected peptide list.')
        self.new_pep_btn = tk.Button(pep_bottom_frame, text="New", command=self.new_pep_prof)
        self.new_pep_btn.grid(row=0,column=3,padx=3)
        new_pep_tip = Hovertip(self.new_pep_btn,'Write a new peptide list.')
        
    def ve_pep_change(self,*args):
        self.pep_entry.delete('1.0',tk.END)
        self.pep_entry.insert(tk.END, self.peptide_profile_dict[self.chosen_pep.get()])
    def del_cur_pep(self):
        if self.chosen_pep.get() == 'None':
            return
        del self.peptide_profile_dict[self.chosen_pep.get()]
        with open('user_data/peptide_profiles.txt','w') as f:
            print(self.peptide_profile_dict,file=f)
        self.pepKeys = list(self.peptide_profile_dict.keys())
        self.pep_entry.delete('1.0',tk.END)
        self.pep_options.destroy()
        self.pep_options = ttk.OptionMenu(self,self.chosen_pep, 'None', *self.pepKeys, command = self.ve_pep_change)
        self.pep_options.grid(row=2,column=1)   
    def edit_cur_pep(self):
        if self.chosen_pep.get() == 'None':
            return
        self.peptide_profile_dict[self.chosen_pep.get()] = self.pep_entry.get('1.0',tk.END).upper().replace('\n','').replace(' ','')
        with open('user_data/peptide_profiles.txt','w') as f:
            print(self.peptide_profile_dict,file=f)
        self.pepKeys = list(self.peptide_profile_dict.keys())
        self.pep_options.destroy()
        self.pep_options = ttk.OptionMenu(self,self.chosen_pep, self.chosen_pep.get(), *self.pepKeys, command = self.ve_pep_change)
        self.pep_options.grid(row=2,column=1)
        




        
    def new_pep_prof(self):
        self.destroy_all()
        self.npep_title_lbl = tk.Label(self, text="New Peptides",font=("Helvetica",10,"bold"))
        self.npep_title_lbl.grid(row=0,column=0,columnspan=2)
        self.return_home_btn = tk.Button(self, text="\u2BAA"+" Back", command=self.editPeptidesProfiles)
        self.return_home_btn.grid(row=0,column=0,sticky='w')
        self.columnconfigure(0,weight=1)
        self.columnconfigure(1,weight=4)
        self.npepname_lbl = tk.Label(self, text="List Name:")
        self.npepname_lbl.grid(row=1,column=0)
        self.npepname_entry = tk.Entry(self,width=23)
        self.npepname_entry.grid(row=1,column=1)
        self.newpep_label = tk.Label(self, text="Protein\nSequence:")
        self.newpep_label.grid(row=2,column=0)
        newpepname_tip = Hovertip(self.npepname_lbl,'Enter the name for the new\npeptide list.')
        self.newpep_entry = st.ScrolledText(self,width=23,height=5,font=("Courier",9))
        self.newpep_entry.grid(row=2,column=1)
        newpep_tip = Hovertip(self.newpep_label,'Type or paste in a new peptide list,\nseparated by commas.')
        self.save_new_pep_btn = tk.Button(self, text="Save", command=self.saveNewPepProf)
        self.save_new_pep_btn.grid(row=3,column=1)
    def saveNewPepProf(self):
        self.peptide_profile_dict[self.npepname_entry.get().replace(' ','_')] = self.newpep_entry.get('1.0',tk.END).upper().replace('\n','').replace(' ','')
        self.pepKeys = list(self.peptide_profile_dict.keys())
        with open('user_data/peptide_profiles.txt','w') as f:
            print(self.peptide_profile_dict,file=f)
        self.save_new_pep_btn.configure(text='Saved',state=tk.DISABLED)
        
        #################
        ##################################
        #######################################################
        
    def proteinTranslate(self):
        self.title("TPM: Translation")
        self.geometry('255x320')
        self.destroy_all()

        #read protein profiles
        self.protein_profile_dict = {}
        self.protein_profile_dict['None'] = ''
        with open('user_data/protein_profiles.txt') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                self.protein_profile_dict[record.id] = record.seq
        self.protKeys = list(self.protein_profile_dict.keys())

        #prompt user parameters
        self.translation_title = tk.Label(self, text="Enter Parameters",font=("Helvetica",10,"bold"))
        self.translation_title.grid(row=0,column=0,columnspan=2)

        self.return_home_btn = tk.Button(self, text="\u2BAA"+" Back", command=self.create_widgets)
        self.return_home_btn.grid(row=0,column=0,sticky='w')
        
        self.columnconfigure(0,weight=1)
        self.columnconfigure(1,weight=4)
        
        self.sr_start_label = tk.Label(self, text="Search Range Start:")
        self.sr_start_label.grid(row=1,column=0)
        self.search_range_start = tk.Entry(self,width=15)
        self.search_range_start.grid(row=1,column=1)
        sr_start_tip = Hovertip(self.sr_start_label,'Position of first nucleotide in\nsearch range.\n(Leave at 0 to search entire genome)')
                                
        self.sr_end_label = tk.Label(self, text="Search Range End:")
        self.sr_end_label.grid(row=2,column=0)
        self.search_range_end = tk.Entry(self,width=15)
        self.search_range_end.grid(row=2,column=1)
        sr_end_tip = Hovertip(self.sr_end_label,'Position of last nucleotide in\nsearch range.\n(Leave at 0 to search entire genome)')
        
        self.mti_label = tk.Label(self, text="Max # Indels:")
        self.mti_label.grid(row=3,column=0)
        self.mti_entry = tk.Entry(self,width=15)
        self.mti_entry.grid(row=3,column=1)
        mti_tip = Hovertip(self.mti_label,'Maximum number of tolderated insertions\nor deletions allowed in translations.\nMust be greater than 0.')

        self.match_label = tk.Label(self, text="Match Length:")
        self.match_label.grid(row=4,column=0)
        self.match_entry = tk.Entry(self,width=15)
        self.match_entry.grid(row=4,column=1)
        match_tip = Hovertip(self.match_label,'Length of terminal residues used for sequence searching.\nTypically set at at least 20 and at most\nhalf the length of the reference protein.')

        self.identity_label = tk.Label(self, text="Identity Filter:")
        self.identity_label.grid(row=5,column=0)
        self.identity_entry = tk.Entry(self,width=15)
        self.identity_entry.grid(row=5,column=1)
        identity_tip = Hovertip(self.identity_label,'Minimum identity of translations relative\nto reference. Set between 0 and 1.')
        
        self.ref_label = tk.Label(self, text="Reference\nSequence:")
        self.ref_label.grid(row=6,column=0)
        self.ref_entry = st.ScrolledText(self,width=20,height=5,font=("Courier",9))
        self.ref_entry.grid(row=6,column=1)
        ref_tip = Hovertip(self.ref_label,'The amino acid sequence to be used\nas a reference for translation. Type\nthe sequence or select a predefined sequence\nfrom the dropdown menu.')
        
        self.chosen_ref = tk.StringVar(self)
        self.ref_options = ttk.OptionMenu(self,self.chosen_ref, 'None', *self.protKeys, command = self.prof_change)
        self.ref_options.grid(row=7,column=1)

        auto_prot_path = ''
        for i in range(len(self.selected_file)-1,0,-1):
            if self.selected_file[i] == '.':
                for j in range(i-1,0,-1):
                    if self.selected_file[j] == '/':
                        break
                    auto_prot_path = self.selected_file[j] + auto_prot_path
                break
        auto_prot_path += '_TranslationReport'        
        self.output_label = tk.Label(self, text="Output File:")
        self.output_label.grid(row=8,column=0)
        self.output_entry = st.ScrolledText(self,width=20,height=0.1,font=("Courier",9))
        self.output_entry.insert(tk.END,auto_prot_path)
        self.output_entry.grid(row=8,column=1)
        output_tip = Hovertip(self.output_label,'Name of the output file, to be saved in\nthe \'Outputs\' folder.')
        self.execute_translation_btn = tk.Button(self, text="Execute", command=self.protExec)
        self.execute_translation_btn.grid(row=9,column=0,columnspan=2,pady=3)

        
    def prof_change(self,*args):
        self.ref_entry.delete('1.0',tk.END)
        self.ref_entry.insert(tk.END, self.protein_profile_dict[self.chosen_ref.get()])

    def protExec(self):
        self.translation_params = []
        try:
            reference = self.ref_entry.get('1.0',tk.END).upper().replace('\n','').replace(' ','')
            search_range = (int(self.search_range_start.get()),int(self.search_range_end.get()))
            match_length = int(self.match_entry.get())
            max_indels = int(self.mti_entry.get())
            output = self.output_entry.get('1.0',tk.END).replace('\n','')
            identity = float(self.identity_entry.get())
            if not (reference.isalpha() and (search_range[0]<search_range[1] or (search_range[0]==0 and search_range[1]==0))
                    and match_length>0 and max_indels>=0 and identity <=1):
                raise Exception('Invalid parameters')
            self.destroy_all()

            self.translation_params = (reference,search_range,match_length,output,
                                           max_indels,identity,self.selected_file)
                
            self.title("TPM: Processing")
            self.geometry('255x255')
            self.progress = ttk.Progressbar(self, orient="horizontal", length=200, mode="indeterminate")
            self.progress.place(relx=0.5,rely=0.4,anchor='center')
            self.progress.start()
            self.process_lbl = tk.Label(self, text="Performing translation, please\nleave this window open...")
            self.process_lbl.place(relx=0.5,rely=0.25,anchor='center')
            translation_thread = threading.Thread(target=self.getProtSeqs)
            translation_thread.daemon = True
            translation_thread.start()
        except:
            self.error_reset()


    def peptideReport(self):
        self.title("TPM: Pep. Rep.")
        self.geometry('255x258')
        self.destroy_all()

        #read peptide profiles from saved text file to peptide_profile_dict
        with open('user_data/peptide_profiles.txt','r') as file:
            self.peptide_profile_dict = ast.literal_eval(file.read())
        self.peptide_profile_dict['None'] = ''
        self.pepKeys = list(self.peptide_profile_dict.keys())

        #take user parameters
        self.report_title = tk.Label(self, text="Enter Parameters",font=("Helvetica",10,"bold"))
        self.report_title.grid(row=0,column=0,columnspan=2)

        self.return_home_btn = tk.Button(self, text="\u2BAA"+" Back", command=self.create_widgets)
        self.return_home_btn.grid(row=0,column=0,sticky='w')
        
        self.columnconfigure(0,weight=1)
        self.columnconfigure(1,weight=4)

        self.enzyme_label = tk.Label(self, text="Enzyme for Digest:")
        self.enzyme_label.grid(row=1,column=0)
        enzymes = ['Trypsin','a-Lytic']
        self.chosen_enz = tk.StringVar(self)
        self.enzyme_choice = ttk.OptionMenu(self,self.chosen_enz, enzymes[0], *enzymes)
        self.enzyme_choice.grid(row=1,column=1)
        enzyme_tip = Hovertip(self.enzyme_label,'Select enzyme to define cleavage\nsites and delineate peptides.')
        
        self.threshold_label = tk.Label(self, text="Min. Abundance:")
        self.threshold_label.grid(row=2,column=0)
        self.threshold_entry = tk.Entry(self,width=15)
        self.threshold_entry.grid(row=2,column=1)
        threshold_tip = Hovertip(self.threshold_label,'Only peptides found at or above this\nvalue will be included in the report.\nMust be between 0-100. You can also\nset to -1 to display the 5 most\nabundant alternate peptides instead.')

        self.pepref_label = tk.Label(self, text="Peptides\nof Interest:")
        self.pepref_label.grid(row=3,column=0)
        self.pepref_entry = st.ScrolledText(self,width=20,height=5,font=("Courier",9))
        self.pepref_entry.grid(row=3,column=1)
        pepref_tip = Hovertip(self.pepref_label,'The peptides that will be searched for\nin the given file. Type them in, separated by\ncommas, or select a list from the dropdown menu.')
        self.chosen_pepref = tk.StringVar(self)
        self.pepref_options = ttk.OptionMenu(self,self.chosen_pepref, 'None', *self.pepKeys, command = self.pep_prof_change)
        self.pepref_options.grid(row=4,column=1)
        
        auto_pep_path = ''
        for i in range(len(self.selected_file)-1,0,-1):
            if self.selected_file[i] == '.':
                for j in range(i-1,0,-1):
                    if self.selected_file[j] == '/':
                        break
                    auto_pep_path = self.selected_file[j] + auto_pep_path
                break
        auto_pep_path += '_PeptideReport'        
        self.output_label = tk.Label(self, text="Output File:")
        self.output_label.grid(row=8,column=0)
        self.output_entry = st.ScrolledText(self,width=20,height=0.1,font=("Courier",9))
        self.output_entry.insert(tk.END,auto_pep_path)
        self.output_entry.grid(row=8,column=1)
        output_tip = Hovertip(self.output_label,'Name of the output file, to be saved in\nthe \'Outputs\' folder.')

        self.pep_execute_translation_btn = tk.Button(self, text="Execute", command=self.pepExec)
        self.pep_execute_translation_btn.grid(row=9,column=0,columnspan=2,pady=3)
        
    def pep_prof_change(self,*args):
        self.pepref_entry.delete('1.0',tk.END)
        self.pepref_entry.insert(tk.END, self.peptide_profile_dict[self.chosen_pepref.get()])

    def pepExec(self):
        self.report_params = []
        try:
            peptides = self.pepref_entry.get('1.0',tk.END).upper().replace('\n','').replace(' ','').split(',')
            min_threshold = int(self.threshold_entry.get())
            digest_enzyme = (self.chosen_enz.get())
            output = self.output_entry.get('1.0',tk.END).replace('\n','')
                
            if not (min_threshold<=100 and len(output)>0 and len(peptides)>0):
                raise Exception('Invalid parameters')
            for pep in peptides:
                if not (pep.isalpha()):
                    raise Exception('Invalid parameters')
            self.destroy_all()
        
            self.report_params = (peptides,digest_enzyme,output,min_threshold)
            self.title("TPM: Processing")
            self.geometry('255x255')
            self.progress = ttk.Progressbar(self, orient="horizontal", length=200, mode="indeterminate")
            self.progress.place(relx=0.5,rely=0.4,anchor='center')
            self.progress.start()
            self.process_lbl = tk.Label(self, text="Performing peptide analysis,\nplease leave this window open...")
            self.process_lbl.place(relx=0.5,rely=0.25,anchor='center')
            report_thread = threading.Thread(target=self.getPepReport)
            report_thread.daemon = True
            report_thread.start()
            
        except:
            self.error_reset()

    
    def getProtSeqs(self):
        paramsList = self.translation_params
        output_path = 'outputs/'+paramsList[3]+'.txt'
        search_range = paramsList[1]
        match_length = paramsList[2]
        reference = paramsList[0][:match_length]
        reference2 = paramsList[0][-match_length:]
        max_indel_tolerance = paramsList[4]
        pct_identity = paramsList[5]
        fasta_path = paramsList[6]

        codon_dict = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W', }

        #pam100,sc0.4
        pam_100_matrix = {'A': {'A': 3, 'R': -3, 'N': -1, 'D': -1, 'C': -3, 'Q': -1, 'E': 0, 'G': 0, 'H': -3, 'I': -1, 'L': -3, 'K': -3, 'M': -2, 'F': -4, 'P': 0, 'S': 1, 'T': 1, 'W': -6, 'Y': -4, 'V': 0}, 'R': {'A': -3, 'R': 6, 'N': -1, 'D': -3, 'C': -4, 'Q': 0, 'E': -3, 'G': -4, 'H': 1, 'I': -2, 'L': -4, 'K': 2, 'M': -1, 'F': -5, 'P': -1, 'S': -1, 'T': -2, 'W': 1, 'Y': -5, 'V': -3}, 'N': {'A': -1, 'R': -1, 'N': 4, 'D': 2, 'C': -5, 'Q': 0, 'E': 0, 'G': -1, 'H': 1, 'I': -2, 'L': -3, 'K': 1, 'M': -3, 'F': -4, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3}, 'D': {'A': -1, 'R': -3, 'N': 2, 'D': 5, 'C': -6, 'Q': 0, 'E': 3, 'G': -1, 'H': -1, 'I': -3, 'L': -5, 'K': -1, 'M': -4, 'F': -7, 'P': -3, 'S': -1, 'T': -1, 'W': -8, 'Y': -5, 'V': -3}, 'C': {'A': -3, 'R': -4, 'N': -5, 'D': -6, 'C': 8, 'Q': -7, 'E': -7, 'G': -4, 'H': -4, 'I': -3, 'L': -7, 'K': -7, 'M': -6, 'F': -6, 'P': -4, 'S': -1, 'T': -3, 'W': -8, 'Y': -1, 'V': -3}, 'Q': {'A': -1, 'R': 0, 'N': 0, 'D': 0, 'C': -7, 'Q': 5, 'E': 2, 'G': -3, 'H': 2, 'I': -3, 'L': -2, 'K': 0, 'M': -1, 'F': -6, 'P': 0, 'S': -2, 'T': -2, 'W': -6, 'Y': -5, 'V': -3}, 'E': {'A': 0, 'R': -3, 'N': 0, 'D': 3, 'C': -7, 'Q': 2, 'E': 5, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -7, 'P': -2, 'S': -1, 'T': -2, 'W': -8, 'Y': -5, 'V': -3}, 'G': {'A': 0, 'R': -4, 'N': -1, 'D': -1, 'C': -4, 'Q': -3, 'E': -1, 'G': 5, 'H': -4, 'I': -4, 'L': -5, 'K': -3, 'M': -4, 'F': -5, 'P': -2, 'S': 0, 'T': -2, 'W': -8, 'Y': -6, 'V': -2}, 'H': {'A': -3, 'R': 1, 'N': 1, 'D': -1, 'C': -4, 'Q': 2, 'E': -1, 'G': -4, 'H': 6, 'I': -4, 'L': -3, 'K': -2, 'M': -4, 'F': -3, 'P': -1, 'S': -2, 'T': -3, 'W': -3, 'Y': -1, 'V': -3}, 'I': {'A': -1, 'R': -2, 'N': -2, 'D': -3, 'C': -3, 'Q': -3, 'E': -3, 'G': -4, 'H': -4, 'I': 5, 'L': 1, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -3, 'T': 0, 'W': -6, 'Y': -2, 'V': 3}, 'L': {'A': -3, 'R': -4, 'N': -3, 'D': -5, 'C': -7, 'Q': -2, 'E': -4, 'G': -5, 'H': -3, 'I': 1, 'L': 5, 'K': -4, 'M': 2, 'F': 0, 'P': -3, 'S': -4, 'T': -3, 'W': -3, 'Y': -3, 'V': 0}, 'K': {'A': -3, 'R': 2, 'N': 1, 'D': -1, 'C': -7, 'Q': 0, 'E': -1, 'G': -3, 'H': -2, 'I': -3, 'L': -4, 'K': 5, 'M': 0, 'F': -6, 'P': -3, 'S': -1, 'T': -1, 'W': -5, 'Y': -5, 'V': -4}, 'M': {'A': -2, 'R': -1, 'N': -3, 'D': -4, 'C': -6, 'Q': -1, 'E': -3, 'G': -4, 'H': -4, 'I': 1, 'L': 2, 'K': 0, 'M': 7, 'F': -1, 'P': -3, 'S': -2, 'T': -1, 'W': -6, 'Y': -4, 'V': 1}, 'F': {'A': -4, 'R': -5, 'N': -4, 'D': -7, 'C': -6, 'Q': -6, 'E': -7, 'G': -5, 'H': -3, 'I': 0, 'L': 0, 'K': -6, 'M': -1, 'F': 7, 'P': -5, 'S': -3, 'T': -4, 'W': -1, 'Y': 4, 'V': -3}, 'P': {'A': 0, 'R': -1, 'N': -2, 'D': -3, 'C': -4, 'Q': 0, 'E': -2, 'G': -2, 'H': -1, 'I': -3, 'L': -3, 'K': -3, 'M': -3, 'F': -5, 'P': 6, 'S': 0, 'T': -1, 'W': -6, 'Y': -6, 'V': -2}, 'S': {'A': 1, 'R': -1, 'N': 1, 'D': -1, 'C': -1, 'Q': -2, 'E': -1, 'G': 0, 'H': -2, 'I': -3, 'L': -4, 'K': -1, 'M': -2, 'F': -3, 'P': 0, 'S': 3, 'T': 1, 'W': -2, 'Y': -3, 'V': -2}, 'T': {'A': 1, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': -2, 'H': -3, 'I': 0, 'L': -3, 'K': -1, 'M': -1, 'F': -4, 'P': -1, 'S': 1, 'T': 4, 'W': -6, 'Y': -3, 'V': 0}, 'W': {'A': -6, 'R': 1, 'N': -4, 'D': -8, 'C': -8, 'Q': -6, 'E': -8, 'G': -8, 'H': -3, 'I': -6, 'L': -3, 'K': -5, 'M': -6, 'F': -1, 'P': -6, 'S': -2, 'T': -6, 'W': 11, 'Y': -2, 'V': -7}, 'Y': {'A': -4, 'R': -5, 'N': -2, 'D': -5, 'C': -1, 'Q': -5, 'E': -5, 'G': -6, 'H': -1, 'I': -2, 'L': -3, 'K': -5, 'M': -4, 'F': 4, 'P': -6, 'S': -3, 'T': -3, 'W': -2, 'Y': 7, 'V': -3}, 'V': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -3, 'Q': -3, 'E': -3, 'G': -2, 'H': -3, 'I': 3, 'L': 0, 'K': -4, 'M': 1, 'F': -3, 'P': -2, 'S': -2, 'T': 0, 'W': -7, 'Y': -3, 'V': 5}}
        #pam50,sc0.5
        pam_50_matrix = {'A': {'A': 4, 'R': -4, 'N': -2, 'D': -1, 'C': -3, 'Q': -2, 'E': -1, 'G': 0, 'H': -4, 'I': -2, 'L': -3, 'K': -4, 'M': -2, 'F': -5, 'P': 0, 'S': 0, 'T': 0, 'W': -7, 'Y': -4, 'V': -1}, 'R': {'A': -4, 'R': 6, 'N': -3, 'D': -5, 'C': -4, 'Q': 0, 'E': -5, 'G': -5, 'H': 0, 'I': -3, 'L': -5, 'K': 1, 'M': -2, 'F': -5, 'P': -2, 'S': -1, 'T': -3, 'W': 0, 'Y': -6, 'V': -4}, 'N': {'A': -2, 'R': -3, 'N': 5, 'D': 2, 'C': -6, 'Q': -1, 'E': -1, 'G': -1, 'H': 1, 'I': -3, 'L': -4, 'K': 0, 'M': -4, 'F': -5, 'P': -3, 'S': 0, 'T': -1, 'W': -5, 'Y': -2, 'V': -4}, 'D': {'A': -1, 'R': -5, 'N': 2, 'D': 5, 'C': -8, 'Q': -1, 'E': 2, 'G': -1, 'H': -2, 'I': -4, 'L': -7, 'K': -2, 'M': -6, 'F': -8, 'P': -4, 'S': -2, 'T': -2, 'W': -8, 'Y': -6, 'V': -4}, 'C': {'A': -3, 'R': -4, 'N': -6, 'D': -8, 'C': 7, 'Q': -8, 'E': -8, 'G': -5, 'H': -4, 'I': -3, 'L': -8, 'K': -8, 'M': -7, 'F': -7, 'P': -4, 'S': -1, 'T': -4, 'W': -9, 'Y': -2, 'V': -3}, 'Q': {'A': -2, 'R': 0, 'N': -1, 'D': -1, 'C': -8, 'Q': 5, 'E': 1, 'G': -4, 'H': 1, 'I': -4, 'L': -3, 'K': -1, 'M': -2, 'F': -7, 'P': -1, 'S': -3, 'T': -3, 'W': -7, 'Y': -6, 'V': -4}, 'E': {'A': -1, 'R': -5, 'N': -1, 'D': 2, 'C': -8, 'Q': 1, 'E': 5, 'G': -2, 'H': -2, 'I': -3, 'L': -5, 'K': -2, 'M': -4, 'F': -8, 'P': -3, 'S': -2, 'T': -3, 'W': -9, 'Y': -5, 'V': -3}, 'G': {'A': 0, 'R': -5, 'N': -1, 'D': -1, 'C': -5, 'Q': -4, 'E': -2, 'G': 4, 'H': -5, 'I': -6, 'L': -6, 'K': -4, 'M': -5, 'F': -5, 'P': -3, 'S': 0, 'T': -3, 'W': -8, 'Y': -7, 'V': -3}, 'H': {'A': -4, 'R': 0, 'N': 1, 'D': -2, 'C': -4, 'Q': 1, 'E': -2, 'G': -5, 'H': 6, 'I': -5, 'L': -3, 'K': -3, 'M': -5, 'F': -3, 'P': -2, 'S': -3, 'T': -4, 'W': -4, 'Y': -1, 'V': -3}, 'I': {'A': -2, 'R': -3, 'N': -3, 'D': -4, 'C': -3, 'Q': -4, 'E': -3, 'G': -6, 'H': -5, 'I': 5, 'L': 0, 'K': -3, 'M': 0, 'F': -1, 'P': -5, 'S': -3, 'T': -1, 'W': -7, 'Y': -3, 'V': 2}, 'L': {'A': -3, 'R': -5, 'N': -4, 'D': -7, 'C': -8, 'Q': -3, 'E': -5, 'G': -6, 'H': -3, 'I': 0, 'L': 4, 'K': -4, 'M': 1, 'F': -1, 'P': -4, 'S': -5, 'T': -4, 'W': -3, 'Y': -4, 'V': -1}, 'K': {'A': -4, 'R': 1, 'N': 0, 'D': -2, 'C': -8, 'Q': -1, 'E': -2, 'G': -4, 'H': -3, 'I': -3, 'L': -4, 'K': 4, 'M': 0, 'F': -8, 'P': -3, 'S': -2, 'T': -1, 'W': -6, 'Y': -5, 'V': -5}, 'M': {'A': -2, 'R': -2, 'N': -4, 'D': -6, 'C': -7, 'Q': -2, 'E': -4, 'G': -5, 'H': -5, 'I': 0, 'L': 1, 'K': 0, 'M': 7, 'F': -2, 'P': -4, 'S': -3, 'T': -2, 'W': -7, 'Y': -6, 'V': 0}, 'F': {'A': -5, 'R': -5, 'N': -5, 'D': -8, 'C': -7, 'Q': -7, 'E': -8, 'G': -5, 'H': -3, 'I': -1, 'L': -1, 'K': -8, 'M': -2, 'F': 6, 'P': -6, 'S': -4, 'T': -5, 'W': -2, 'Y': 2, 'V': -4}, 'P': {'A': 0, 'R': -2, 'N': -3, 'D': -4, 'C': -4, 'Q': -1, 'E': -3, 'G': -3, 'H': -2, 'I': -5, 'L': -4, 'K': -3, 'M': -4, 'F': -6, 'P': 5, 'S': 0, 'T': -2, 'W': -8, 'Y': -7, 'V': -3}, 'S': {'A': 0, 'R': -1, 'N': 0, 'D': -2, 'C': -1, 'Q': -3, 'E': -2, 'G': 0, 'H': -3, 'I': -3, 'L': -5, 'K': -2, 'M': -3, 'F': -4, 'P': 0, 'S': 4, 'T': 1, 'W': -3, 'Y': -4, 'V': -3}, 'T': {'A': 0, 'R': -3, 'N': -1, 'D': -2, 'C': -4, 'Q': -3, 'E': -3, 'G': -3, 'H': -4, 'I': -1, 'L': -4, 'K': -1, 'M': -2, 'F': -5, 'P': -2, 'S': 1, 'T': 4, 'W': -7, 'Y': -3, 'V': -1}, 'W': {'A': -7, 'R': 0, 'N': -5, 'D': -8, 'C': -9, 'Q': -7, 'E': -9, 'G': -8, 'H': -4, 'I': -7, 'L': -3, 'K': -6, 'M': -7, 'F': -2, 'P': -8, 'S': -3, 'T': -7, 'W': 9, 'Y': -3, 'V': -9}, 'Y': {'A': -4, 'R': -6, 'N': -2, 'D': -6, 'C': -2, 'Q': -6, 'E': -5, 'G': -7, 'H': -1, 'I': -3, 'L': -4, 'K': -5, 'M': -6, 'F': 2, 'P': -7, 'S': -4, 'T': -3, 'W': -3, 'Y': 6, 'V': -4}, 'V': {'A': -1, 'R': -4, 'N': -4, 'D': -4, 'C': -3, 'Q': -4, 'E': -3, 'G': -3, 'H': -3, 'I': 2, 'L': -1, 'K': -5, 'M': 0, 'F': -4, 'P': -3, 'S': -3, 'T': -1, 'W': -9, 'Y': -4, 'V': 5}}

        def limitChars(seq,limiter):
            output = ''
            for i in range(len(seq)):
                if i%limiter == 0 and i != 0:
                    output += '\n'
                output += seq[i]
            return output
    
        #find all proteins within genome (i.e. methionine to stop codon) at specified register
        def findProteins(genome, offset):
            prot = ''
            for i in range(offset, len(genome), 3):
                try:
                    prot += codon_dict[genome[i:i + 3]]
                except KeyError:
                    prot += 'X'
            prots_list = []
            for i in range(len(prot)):
                if prot[i] == 'M':
                    thisone = ''
                    for j in range(i, len(prot)):
                        if prot[j] == '_':
                            if len(thisone) <= max_indel_tolerance+len(paramsList[0]) and len(thisone) >= len(paramsList[0])-max_indel_tolerance:
                                prots_list.append(thisone)
                            break
                        thisone += prot[j]
            return prots_list

        def similarity(s1, s2):
            rows, cols = (len(s1) + 1, len(s2) + 1)
            dp = [[(0,False) for i in range(cols)] for j in range(rows)]

            #affine gap penalties
            gap_open = -12
            gap_extend = -3
        
            for i in range(len(dp)):
                for j in range(len(dp[0])):
                    if i == 0:
                        dp[i][j] = (-1*j,False)
                    if j == 0:
                        dp[i][j] = (-1*i,False)
                    if (i > 0 and j > 0):
                        try:
                            match_score = pam_50_matrix[s1[i-1]][s2[j-1]]
                        except:
                            #handles unresolves bases/nucleotides: how forgiving should we be?
                            match_score = -3
                        '''
                        #constant match/mismatch
                        match_score = -3
                        if s1[i - 1] == s2[j - 1]:
                            match_score = 3
                        '''
                        match_res = dp[i-1][j-1][0] + match_score
                        up_gap = dp[i-1][j][0]+gap_extend
                        if not dp[i-1][j][1]:
                            up_gap = dp[i-1][j][0]+gap_open
                        left_gap = dp[i][j-1][0]+gap_extend
                        if not dp[i-1][j][1]:
                            left_gap = dp[i][j-1][0]+gap_open
                        if (match_res >= up_gap) and (match_res >= left_gap):
                            dp[i][j] = (match_res,False)
                        elif (up_gap >= left_gap):
                            dp[i][j] = (up_gap,True)
                        else:
                            dp[i][j] = (left_gap,False)
            '''
            #Display grid
            def lfy(A):
                if A:
                    return 'T'
                return 'F'
                    
            for line in dp:
                for item in line:
                    if str(item[0])[0] != '-':
                        print('(+'+str(item[0])+', '+lfy(item[1])+')',end = '')
                    else:
                        print('('+str(item[0])+', '+lfy(item[1])+')',end = '')
                print('\n')
            '''
            return dp[len(dp) - 1][len(dp[0]) - 1][0]

        def expected_score(ref,itr,sim):
            self_match = similarity(ref,ref)
            aa_list = 'ARNDCQEGHILKMFPSTWYV'
            sims = []
            for i in range(itr):
                prot_rand = ''
                for i in range(len(ref)):
                    ind = randint(0,19)
                    prot_rand += aa_list[ind]
                sims.append(similarity(prot_rand,ref))
            return (sim)*(self_match) + (1-sim)*(sum(sims)/itr)
        
        avg_exp_score_1,avg_exp_score_2 = expected_score(reference,250,pct_identity),expected_score(reference2,250,pct_identity)
      
        def findTarget(genome,reference):
            maxProt = (avg_exp_score_1,avg_exp_score_2, '')
            #maxProt = (-1,-1,'')
            for i in range(3):
                these_prots = findProteins(genome, i)
                for prot in these_prots:
                    score_n = similarity(prot[0:len(reference)], reference)
                    score_c = similarity(prot[-len(reference2):],reference2)
                    if (score_n>=maxProt[0]) and (score_c>=maxProt[1]):
                        maxProt = (score_n,score_c, prot)
            return maxProt

        startTime = time.time()
        bad_ones = []
        numSeqs = 0
        full_thing = ''
        #write final output here
        with open(output_path, 'w') as f:
            global today
            print('Genome to protein output for file:\n' + fasta_path + '.\nGenerated on '+ str(today) + '.\n',file=f)
            print('Parameters:\nSearch range = ' + str(search_range) +'\nMatch length = ' + str(match_length) +'\nMax. tolerated indels = '+str(max_indel_tolerance)+'\nIdentity threshold = '+str(pct_identity)+'\nReference:\n'+limitChars(str(paramsList[0]),90)+'\n',file=f)

            with open(fasta_path) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    numSeqs += 1
                    if search_range == (0,0):
                        bestProtein = findTarget(record.seq.upper(), reference)
                    else:
                        try:
                            bestProtein = findTarget(record.seq[search_range[0]:search_range[1]].upper(), reference)
                        except:
                            bestProtein = findTarget(record.seq[search_range[0]:len(record.seq)].upper(), reference)
                    if len(bestProtein[2]) < (len(paramsList[0])*0.97) or len(bestProtein[2]) > (len(paramsList[0])*1.03):
                        bad_ones.append(record.id)
                    else:
                        '''
                        print('>'+record.id, file=f)
                        print(bestProtein[1], file=f)
                        '''
                        full_thing += '\n>'+record.id+'\n'+limitChars(bestProtein[2],100)
            print(str(len(bad_ones))+' sequence(s) where target protein could not be found (potentially obscured start/stop codons):\n'+str(bad_ones),file=f)

            endTime = time.time()
            total = round(endTime-startTime,2)
            print('\nExecution time: '+str(timedelta(seconds=int(total)))+' for '+str(numSeqs)+' genome(s).\n\n*********************************************************',file=f)
            print('\n'+full_thing,file=f)
            
        self.destroy_all()
        self.title('TPM: Complete')
        self.geometry('255x185')
        self.finished = tk.Label(self, text="Execution complete!\nOutput file saved as:")
        self.finished.pack(pady=3)
        self.finished = tk.Label(self, text=limitChars(output_path,30),font=("Courier",9))
        self.finished.pack(pady=3)
        self.finished_btn = tk.Button(self, text="Return Home", command=self.create_widgets)
        self.finished_btn.pack(pady=3)


    def getPepReport(self):
        protParamz = self.report_params
        output_path = 'outputs/'+protParamz[2]+'.txt'
        enz = protParamz[1]
        pepsList = protParamz[0]
        dist_algo = 'NW'
        num_freqs = 5
        pctExclude = float(protParamz[3])
        fasta_path = self.selected_file
        
        if pctExclude <= 0:
            pctExclude = -1
        

        min_pep_length, max_pep_length = 10000, 0
        for i in range(len(pepsList)):
            if len(pepsList[i]) >= max_pep_length:
                max_pep_length = len(pepsList[i])
            if len(pepsList[i]) <= min_pep_length:
                min_pep_length = len(pepsList[i])

        min_pep_length -= 3
        max_pep_length += 3

        def limitChars(seq,limiter):
            output = ''
            for i in range(len(seq)):
                if i%limiter == 0 and i != 0:
                    output += '\n'
                output += seq[i]
            return output

        # ***********enzyme functions********************
        def tryptic(seq, minLn=min_pep_length, maxLn=max_pep_length):
            peps = []
            for i in range(len(seq)):
                if seq[i] == 'K' or seq[i] == 'R':
                    thisPep = ''
                    for j in range(i + 1, len(seq)):
                        thisPep += seq[j]
                        if seq[j] == 'K' or seq[j] == 'R':
                            i += j
                            if (len(thisPep) >= minLn and len(thisPep) <= maxLn):
                                peps.append(thisPep)
                            break
            return peps

        def alp(seq, minLn=min_pep_length, maxLn=max_pep_length):
            peps = []
            for i in range(len(seq)):
                if seq[i] == 'T' or seq[i] == 'A' or seq[i] == 'S' or seq[i] == 'V':
                    thisPep = ''
                    for j in range(i + 1, len(seq)):
                        thisPep += seq[j]
                        if seq[j] == 'T' or seq[j] == 'A' or seq[j] == 'S' or seq[j] == 'V':
                            i += j
                            if (len(thisPep) >= minLn and len(thisPep) <= maxLn):
                                peps.append(thisPep)
                            break
            return peps

        # *************distance functions*************
        def hamming(s1, s2):
            ham = abs(len(s1) - len(s2))
            if len(s1) > len(s2):
                for i in range(len(s2)):
                    if s1[i] != s2[i]:
                        ham += 1
                return ham
            for i in range(len(s1)):
                if s1[i] != s2[i]:
                    ham += 1
            return ham

        def nw(s1, s2):
            rows, cols = (len(s1) + 1, len(s2) + 1)
            dp = [[0 for i in range(cols)] for j in range(rows)]
            for i in range(len(dp)):
                for j in range(len(dp[0])):
                    if i == 0:
                        dp[i][j] = -1 * j
                    if j == 0:
                        dp[i][j] = -1 * i
                    if (i > 0 and j > 0):
                        ms = -1
                        if s1[i - 1] == s2[j - 1]:
                            ms = 1
                        dp[i][j] = max((dp[i - 1][j - 1] + ms), (dp[i - 1][j] - 1), (dp[i][j - 1] - 1))
            return len(s2) - dp[len(dp) - 1][len(dp[0]) - 1]

        # ***************checking peptides************************
        def checking(peps, alg, interest):
            presence = [1] * len(interest)
            for i in range(len(interest)):
                if interest[i] not in peps:
                    min_dist = len(interest[i]) + 1
                    optimal_alt = ''
                    for pep in peps:
                        if alg == 'Hamming':
                            this_dist = hamming(interest[i], pep)
                        if alg == 'NW':
                            this_dist = nw(interest[i], pep)
                        if this_dist <= min_dist:
                            min_dist = this_dist
                            optimal_alt = pep
                    presence[i] = optimal_alt
            return presence

        # ************************driver************************
        startTime = time.time()

        alts = [{} for i in range(1) for j in range(len(pepsList))]
        pep_freq_tally = [0] * len(pepsList)
        details = ''
        num_seqs = 0
        with open(fasta_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                num_seqs += 1
                if enz == 'Trypsin':
                    pepper = tryptic(record.seq)
                if enz == 'A-LyticPrts':
                    pepper = alp(record.seq)
                thisList = checking(pepper, dist_algo, pepsList)
                details += '>' + record.id + '\n' + str(thisList) + '\n\n'
                for i in range(len(thisList)):
                    if thisList[i] == 1:
                        pep_freq_tally[i] += 1
                    else:
                        if thisList[i] not in alts[i]:
                            alts[i][thisList[i]] = 1
                        else:
                            alts[i][thisList[i]] += 1

        def sortAlts(altsList, num):
            answer = []
            for i in range(num_freqs):
                if len(altsList) == 0:
                    break
                maxFreq = ('', 0)
                for i in range(len(altsList)):
                    if altsList[i][1] >= maxFreq[1]:
                        maxFreq = altsList[i]
                        ko = i
                del altsList[ko]
                answer.append(maxFreq)
            return answer

        def sortAltsPct(altsList,pct):
            answer = []
            for i in range(len(altsList)):
                if (altsList[i][1]/num_seqs)*100 >= pct:
                    answer.append(altsList[i])
            return answer

        endTime = time.time()
        total_time = round(endTime - startTime,2)

        with open(output_path, 'w') as f:
            print('Peptide report for file: ' + fasta_path + ', generated on ' + today + '.', file=f)
            print('\nExecution time: ' + str(timedelta(seconds=int(total_time))) + ' for ' + str(num_seqs) + ' sequence(s).\n', file=f)
            for i in range(len(pepsList)):
                print(pepsList[i] + ': ' + str(pep_freq_tally[i]) + ', ' + str(100 * pep_freq_tally[i] / num_seqs)[
                                                                           0:5] + '%\n', file=f)
            if pctExclude == -1:
                print('\n' + str(num_freqs) + ' most frequent alternate peptides found in sequences:', file=f)
            else:
                print('Most frequent alternate peptides found in sequences at >= '+str(protParamz[3])+'% abundance:', file=f)
            for i in range(len(pepsList)):
                if len(alts[i]) != 0:
                    print(pepsList[i], file=f)
                    if pctExclude == -1:
                        these_alts = []
                        for pep in alts[i]:
                            these_alts.append((pep, alts[i][pep]))
                        print(str(sortAlts(these_alts, num_freqs)) + '\n', file=f)
                    else:
                        these_alts = []
                        for pep in alts[i]:
                            these_alts.append((pep, alts[i][pep]))
                        print(str(sortAltsPct(these_alts, pctExclude)) + '\n', file=f)

            print('\n*************************************************\n', file=f)
            print(details, file=f)
        self.destroy_all()
        self.title('TPM: Complete')
        self.geometry('255x185')
        self.finished = tk.Label(self, text="Execution complete!\nOutput file saved as:")
        self.finished.pack(pady=3)
        self.finished = tk.Label(self, text=limitChars(output_path,30),font=("Courier",9))
        self.finished.pack(pady=3)
        self.finished_btn = tk.Button(self, text="Return Home", command=self.create_widgets)
        self.finished_btn.pack(pady=3)


if __name__ == "__main__":
    app = App()
    app.mainloop()
