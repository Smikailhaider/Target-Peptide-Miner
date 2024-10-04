'''
updated file containing newer, working methods/gui

def prot_to_report (line 730):
-takes in fasta file containing amino acid sequences, creates output file containing peptide report

def getProtSeqs (line 627)
-takes in genomes, outputs protein sequences based on reference sequence

def genome_to_pep (line 423)
-takes in genome, creates file containing peptide report

run this, it will create any missing directories used to store files
makes it a little easier for an end user
'''


import os
from tkinter import *
from tkinter import filedialog,scrolledtext
from Bio import SeqIO
import time
from datetime import date
import ast
today = str(date.today())


def validateDirectories():
    return (os.path.isdir('sequences') and os.path.isdir('outputs') and
            os.path.exists('protein_profiles.txt') and os.path.exists('peptide_profiles.txt'))
def buildDirectories():
    for path in ('sequences','outputs'):
        if not os.path.isdir(path):
            os.mkdir(path)
    for file in ('protein_profiles.txt','peptide_profiles.txt'):
        if not os.path.exists(file):
            with open(file,'w') as f:
                pass
    
   
#set up main window
main_display = Tk()
main_display.title('Peptide Finder')
main_display.geometry('255x255')

#create file selector
selected_file = ''
label_1 = Label(main_display, text = '         select the location of the sequences         ')
label_1.grid(row=0,column=0)

def browseFiles():
    for widget in main_display.grid_slaves(row=2):
        widget.configure(text='')
    global selected_file
    selected_file = filedialog.askopenfilename(initialdir = os.getcwd()+'/sequences',
                                          title = "Select a File",
                                          filetypes = (("FASTA files",
                                                        "*.fasta*"),("Text files",
                                                        "*.txt*")
                                                       ))
    display_file = ''
    per_line = 35
    for i in range(0,len(selected_file),per_line):
        if ((i+per_line) > len(selected_file)):
            display_file += selected_file[i:]
        else:
            display_file += selected_file[i:i+per_line]+'\n'
    conf_lbl = Label(main_display, text = 'selected file: \n'+display_file)
    conf_lbl.configure(font = ("Courier", 9, "bold"))
    conf_lbl.grid(column=0,row=2)
    if len(selected_file)>2:
        genome_to_pep_btn['state'] = 'active'
        genome_to_prot_btn['state'] = 'active'
        prot_to_pep_btn['state'] = 'active'
button_explore = Button(main_display, text = "Browse Files",command = browseFiles)
button_explore.grid(row=1,column=0)

conf_lbl = Label(main_display, text = 'no file selected yet\n')
conf_lbl.grid(column=0,row=2)

options_lbl = Label(main_display, text = 'select desired operation')
options_lbl.grid(row=4,column=0)

#method to gray out all buttons (used when execution is in progress)
def grayAll():
    button_explore['state'] = 'disabled'
    genome_to_pep_btn['state'] = 'disabled'
    genome_to_prot_btn['state'] = 'disabled'
    prot_to_pep_btn['state'] = 'disabled'
#show error message, clear window, restart
def errorReset(win):
    clear_frame(win)
    final_lbl = Label(win, text='''An error was encountered. Please ensure
that you are using the correct sequences/files''',fg='red')
    final_lbl.grid(row=8,column=0)
def clear_frame(frame):
    for widget in frame.winfo_children():
        widget.destroy()
def hide_all_widgets(root):
    for widget in root.winfo_children():
        widget.grid_forget()
def show_all_widgets(root):
    for widget in root.winfo_children():
        widget.grid()
    
def show_output_file(frame,file):
    clear_frame(frame)
    final_lbl = Label(frame, text='')
    display_file = ''
    per_line = 35
    for i in range(0,len(file),per_line):
        if ((i+per_line) > len(file)):
            display_file += file[i:]
        else:
            display_file += file[i:i+per_line]+'\n'
    final_lbl.config(text='successfully executed. file saved as:\n'+display_file)
    final_lbl.configure(font = ("Courier", 9, "bold"))
    final_lbl.grid(row=0,column=0)
    return
    

#create button/options for genome -> peptide report
getPepParams = []
def genPepClick():

    #***********************************************
    #read peptide profiles from saved text file to peptide_profile_dict
    with open('peptide_profiles.txt','r') as file:
        peptide_profile_dict = ast.literal_eval(file.read())
    #***********************************************
    
    win = Toplevel()
    win.geometry('270x290')
    win.wm_title("Input Parameters")
    
    l = Label(win, text="Enter search parameters")
    l.grid(row=0, column=0)
    
    sr_text = Entry(win)
    sr_text.insert(0,'20001')
    sr_text.grid(row=1,column=0)
    sr_start_lbl = Label(win,text='search range start')
    sr_start_lbl.grid(row=1,column=1)

    sr2_text = Entry(win)
    sr2_text.insert(0,'25000')
    sr2_text.grid(row=2,column=0)
    sr2_lbl = Label(win,text='search range finish')
    sr2_lbl.grid(row=2,column=1)

    tl_lbl = Label(win,text='distance algorithm')
    tl_lbl.grid(row=3,column=1)
    algo_options = ['hamming','levenshtein']
    alg = StringVar(win)
    alg.set(algo_options[0])
    w = OptionMenu(win,alg,*algo_options)
    w.grid(row=3,column=0)

    pl_area = scrolledtext.ScrolledText(win, wrap=WORD, 
                                      width=14, height=4, pady=2,padx=2)
    pl_area.grid(column=0, row=4)
    pl_lbl = Label(win,text='enter peptides of interest\nseparate with commas')
    pl_lbl.grid(row=4,column=1)

    def changePl(*args):
        pl_area.delete('1.0',END)
        pl_area.insert(END, peptide_profile_dict[prof.get()])
    
    pep_prof_user_selected = ''
    pp_lbl = Label(win,text='select user defined\npeptide profile')
    pp_lbl.grid(row=6,column=1)
    pep_prof_options = list(peptide_profile_dict.keys())
    prof = StringVar(win)
    prof.set('none')
    w2 = OptionMenu(win,prof,*pep_prof_options,command=changePl)
    w2.grid(row=6,column=0)

    auto_report_path = ''
    for i in range(len(selected_file)-1,0,-1):
        if selected_file[i] == '.':
            for j in range(i-1,0,-1):
                if selected_file[j] == '/':
                    break
                auto_report_path = selected_file[j] + auto_report_path
            break
    auto_report_path += '_report.txt'
    report_area = scrolledtext.ScrolledText(win, wrap=WORD, 
                                      width=14, height=1.5, pady=2,padx=2)
    report_area.grid(column=0, row=7)
    report_area.insert(END, auto_report_path)
    report_lbl = Label(win,text='what to name the\nreport file')
    report_lbl.grid(row=7,column=1)


    confirm_btn = Button(win, text="Confirm Parameters", command=lambda:[doPepReport(),execReport()])
    confirm_btn.grid(row=8, column=0)
    def doPepReport():
        getPepParams.append(pl_area.get('1.0', END).replace('\n', '').replace(' ','').split(','))
        getPepParams.append((int(sr_text.get()),int(sr2_text.get())))
        getPepParams.append(alg.get())
        getPepParams.append(report_area.get('1.0', END).replace('\n', '').replace(' ',''))
    def execReport():
        try:
            genome_to_pep(selected_file,getPepParams)
            show_output_file(win,getPepParams[3])
        except:
            errorReset(win)
    return
    
genome_to_pep_btn = Button(main_display, text = 'genome -> peptide report' ,fg = 'red', command=genPepClick)
genome_to_pep_btn['state'] = 'disabled'
genome_to_pep_btn.grid(row=5,column=0)





#driver for genome -> protein, popup window to specify method parameters
getProtParams = []
def getProt():
    
    #read protein profiles
    with open('peptide_profiles.txt','r') as file:
        peptide_profile_dict = ast.literal_eval(file.read())
    protein_profile_dict['none'] = ''
    
    win = Toplevel()
    win.geometry('270x290')
    win.wm_title("Input Parameters")
    l = Label(win, text="Enter search parameters")
    l.grid(row=0, column=0)
    b = Button(win, text="Execute", command=lambda: [updateProtParams(),protExec()])
    b.grid(row=11, column=0)
    sr_text = Entry(win)
    sr_text.insert(0,'18000')
    sr_text.grid(row=1,column=0)
    sr_start_lbl = Label(win,text='search range start')
    sr_start_lbl.grid(row=1,column=1)

    sr2_text = Entry(win)
    sr2_text.insert(0,'26500')
    sr2_text.grid(row=2,column=0)
    sr2_lbl = Label(win,text='search range finish')
    sr2_lbl.grid(row=2,column=1)

    tl_text = Entry(win)
    tl_text.insert(0,'75')
    tl_text.grid(row=3,column=0)
    tl_lbl = Label(win,text='match length\n(<=75 recommended)')
    tl_lbl.grid(row=3,column=1)

    rf_area = scrolledtext.ScrolledText(win, wrap=WORD, 
                                      width=14, height=4, pady=2,padx=2)
    rf_area.grid(column=0, row=4)
    rf_area.insert(END,'')
    rf_lbl = Label(win,text='reference sequence')
    rf_lbl.grid(row=4,column=1)

    def changePl(*args):
        rf_area.delete('1.0',END)
        rf_area.insert(END, protein_profile_dict[prof.get()])
    
    prot_prof_user_selected = ''
    pp_lbl = Label(win,text='select user defined\nprotein sequence')
    pp_lbl.grid(row=6,column=1)
    prot_prof_options = list(protein_profile_dict.keys())
    prof = StringVar(win)
    prof.set('none')
    w2 = OptionMenu(win,prof,*prot_prof_options,command=changePl)
    w2.grid(row=6,column=0)

    auto_prot_path = ''
    for i in range(len(selected_file)-1,0,-1):
        if selected_file[i] == '.':
            for j in range(i-1,0,-1):
                if selected_file[j] == '/':
                    break
                auto_prot_path = selected_file[j] + auto_prot_path
            break
    auto_prot_path += '_toProtein.txt'
    report_area = scrolledtext.ScrolledText(win, wrap=WORD, 
                                      width=14, height=1.5, pady=2,padx=2)
    report_area.grid(column=0, row=9)
    report_area.insert(END, auto_prot_path)
    report_lbl = Label(win,text='what to name the\noutput file')
    report_lbl.grid(row=9,column=1)
    def updateProtParams():
        getProtParams.append(rf_area.get('1.0', END).replace('\n', '').replace(' ',''))
        getProtParams.append((int(sr_text.get()),int(sr2_text.get())))
        getProtParams.append(int(tl_text.get()))
        getProtParams.append(report_area.get('1.0', END).replace('\n', '').replace(' ',''))
    def protExec():
        
        try:
            getProtSeqs(selected_file,getProtParams)
            show_output_file(win,getProtParams[3])
        except:
            errorReset(win)

genome_to_prot_btn = Button(main_display, text = 'genome -> target protein' ,fg = 'red', command=getProt)
genome_to_prot_btn['state'] = 'disabled'
genome_to_prot_btn.grid(row=6,column=0)


#driver for protein -> peps, popup window to specify method parameters
getProtPepParams = []
def getProtReport():

    #***********************************************
    #read peptide profiles from saved text file to peptide_profile_dict
    with open('peptide_profiles.txt','r') as file:
        peptide_profile_dict = ast.literal_eval(file.read())
    peptide_profile_dict['none'] = ''
    #***********************************************
    
    win = Toplevel()
    win.geometry('270x290')
    win.wm_title("Input Parameters")
    l = Label(win, text="Enter report parameters")
    l.grid(row=0, column=0)

    tl_lbl = Label(win, text='distance algorithm')
    tl_lbl.grid(row=3, column=1)
    algo_options = ['hamming', 'levenshtein']
    alg = StringVar(win)
    alg.set(algo_options[0])
    w = OptionMenu(win, alg, *algo_options)
    w.grid(row=3, column=0)

    tl_lbl2 = Label(win, text='enzyme')
    tl_lbl2.grid(row=4, column=1)
    enz_options = ['trypsin', 'alp']
    enz = StringVar(win)
    enz.set(enz_options[0])
    w2 = OptionMenu(win, enz, *enz_options)
    w2.grid(row=4, column=0)

    global pep_prof_user_selected
    pl_area = scrolledtext.ScrolledText(win, wrap=WORD,
                                        width=14, height=4, pady=2, padx=2)
    pl_area.grid(column=0, row=5)
    pl_lbl = Label(win, text='enter peptides of\ninterest separate with\ncommas no spaces')
    pl_lbl.grid(row=5, column=1)

    def changePl(*args):
        pl_area.delete('1.0',END)
        pl_area.insert(END, peptide_profile_dict[prof.get()])
    
    pep_prof_user_selected = ''
    pp_lbl = Label(win,text='select user defined\npeptide profile')
    pp_lbl.grid(row=6,column=1)
    pep_prof_options = list(peptide_profile_dict.keys())
    prof = StringVar(win)
    prof.set('none')
    w2 = OptionMenu(win,prof,*pep_prof_options,command=changePl)
    w2.grid(row=6,column=0)

    note_lbl = Label(win, text='*leave a note maybe here')
    note_lbl.grid(row=7, column=0)

    b = Button(win, text="Confirm Parameters",command=lambda: [updateProtParams(),prot_Exec()])
    b.grid(row=10, column=0)

    auto_report_path = ''
    for i in range(len(selected_file)-1,0,-1):
        if selected_file[i] == '.':
            for j in range(i-1,0,-1):
                if selected_file[j] == '/':
                    break
                auto_report_path = selected_file[j] + auto_report_path
            break
    auto_report_path += '_report.txt'
    report_area = scrolledtext.ScrolledText(win, wrap=WORD, 
                                      width=14, height=1.5, pady=2,padx=2)
    report_area.grid(column=0, row=8)
    report_area.insert(END, auto_report_path)
    report_lbl = Label(win,text='what to name the\nreport file')
    report_lbl.grid(row=8,column=1)

    def updateProtParams():
        getProtPepParams.append(pl_area.get('1.0', END).replace('\n', '').replace(' ','').split(','))
        getProtPepParams.append(alg.get())
        getProtPepParams.append(enz.get())
        getProtPepParams.append(report_area.get('1.0', END).replace('\n', '').replace(' ',''))
    def prot_Exec():
        try:
            prot_to_report(selected_file,getProtPepParams)
            show_output_file(win,getProtPepParams[3])
        except:
            errorReset(win)
    return
prot_to_pep_btn = Button(main_display, text = 'protein -> peptide report' ,fg = 'red', command=getProtReport)
prot_to_pep_btn['state'] = 'disabled'
prot_to_pep_btn.grid(row=7,column=0)



#########################################################################################
#########################################################################################
def genome_to_pep(fasta_path,paramz):
    report_path = 'outputs/'+paramz[3]+'.txt'
    #list of peptides of interest
    pepsList = paramz[0]
    #find smallest and biggest peptides
    min_pep_length,max_pep_length = 1000,0
    for i in range(len(pepsList)):
        if len(pepsList[i]) >= max_pep_length:
            max_pep_length = len(pepsList[i]) 
        if len(pepsList[i])<= min_pep_length:
            min_pep_length = len(pepsList[i])
    min_pep_length -= 3
    max_pep_length += 3
    #not necessary, but nice to start at a multiple of 3 
    search_range = paramz[1]
    #enter hamming or levenshtein
    match_algo = paramz[2]

    #maps codons to corresponding amino acids
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

    #translates a nucleotide sequence to aa sequence
    def translate(genome, offset):
        prot = ''
        for i in range(offset, len(genome), 3):
            try:
                prot += codon_dict[genome[i:i + 3]]
            except KeyError:
                #unresolved nucleotide in sequence: just put an X as the next residue
                prot += 'X'
        return prot

    def scoreReg(reg,seq,peps):
        resultant = tryptic(translate(seq,reg))
        score = 0
        presence_list = [0]*len(peps)
        for i in range(len(peps)):
            thisMatch = bestMatch(resultant,peps[i],True)
            if thisMatch == 1:
                presence_list[i] = 1
                score += len(peps[i])
            else:
                presence_list[i] = thisMatch[1]
                score += (len(peps[i])-thisMatch[0])
        return (presence_list,score)            

    #takes a nucelotide sequence (seq) and list of peptide of interest (peps)
    #returns tuple of list containing peptide presence/alternates, and the appropriate reading frame (reg)
    def findPeps(seq, peps):
        maxScore,reg = -1,-1
        best_peps = []
        for i in range(3):
            thisReg = scoreReg(i,seq,peps)
            if thisReg[1] >= maxScore:
                maxScore = thisReg[1]
                best_peps = thisReg[0]
                reg = i
        return (best_peps,reg)

    #take in a aa sequence, cleave at c-term lysine and arginine
    def tryptic(prot):
        peps = []
        for i in range(len(prot)):
            if prot[i] == 'K' or prot[i] == 'R':
                thisOne = ''
                for j in range(i+1,len(prot)):
                    if prot[j] == 'K' or prot[j] == 'R':
                        thisOne += prot[j]
                        if (len(thisOne)>=min_pep_length and len(thisOne)<=max_pep_length and ('_' not in thisOne)):
                            peps.append(thisOne)
                        i = j
                        break
                    thisOne += prot[j]
        return peps

    #hamming distance between two sequences/strings
    def hamming(s1,s2):
        ham = abs(len(s1)-len(s2))
        if len(s1)>len(s2):
            for i in range(len(s2)):
                if s1[i] != s2[i]:
                    ham += 1
            return ham
        for i in range(len(s1)):
                if s1[i] != s2[i]:
                    ham += 1
        return ham

    #levenshtein distance between two strings (i think? based on nw) analog to hamming here
    def levenshtein(s1, s2):
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
        return len(s2)-dp[len(dp) - 1][len(dp[0]) - 1]

    #takes in a tryptic peptide list and one peptide of interest. returns best match with specified algorithm
    def bestMatch(trp,pep,full=False):
        if pep in trp:
            return 1
        min_dist = (len(pep) + 1,'')
        for match in trp:
            if match_algo == 'hamming':
                thisDist = hamming(match,pep)
            if match_algo == 'levenshtein':
                thisDist = levenshtein(match,pep)
            if thisDist <= min_dist[0]:
                min_dist = (thisDist,match)
        if full:
            return min_dist
        return min_dist[1]

    #creates a list of dicts to contain alt peptides + frequencies for each peptide in list
    alts = [{} for i in range(1) for j in range(len(pepsList))]
    startTime = time.time()
    registers = []
    numseqs = 0
    pep_freq_tally = [0] * len(pepsList)
    full_thing = ''
    with open(fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            numseqs += 1
            full_thing += '\n>'+record.id
            if search_range == (0,0):
                pep_summary = findPeps((record.seq).upper(), pepsList)
            else:
                pep_summary = findPeps(record.seq[search_range[0]:search_range[1]].upper(),pepsList)
            full_thing += '\n'+str(pep_summary[0])
            full_thing += '\nregister: '+str(pep_summary[1])+'\n'
            for i in range(len(pepsList)):
                if type(pep_summary[0][i]) is int:
                    pep_freq_tally[i] += pep_summary[0][i]
                else:
                    if pep_summary[0][i] not in alts[i]:
                        alts[i][pep_summary[0][i]] = 1
                    else:
                        alts[i][pep_summary[0][i]] += 1
            registers.append(pep_summary[1])

    num_freqs = 5
    def sortAlts(altsList,num):
        answer = []
        for i in range(num_freqs):
            if len(altsList) == 0:
                break
            maxFreq = ('',0)
            for i in range(len(altsList)):
                if altsList[i][1] >= maxFreq[1]:
                    maxFreq = altsList[i]
                    ko = i
            del altsList[ko]
            answer.append(maxFreq)
        return answer
    endTime = time.time()
    total = str(endTime-startTime)

    
    #writes the actual report, copying data from 'full_thing' str
    with open(report_path,'w') as f:
        global today
        print('genome_to_pepReport report output for file ' + fasta_path + '. \nGenerated on '+ str(today) + '\n',file=f)
        print('search range: '+ str(search_range),file=f)
        print('\nExecution time: ' + total + ' seconds\nnumber of genomes: ' + str(numseqs) + '. peptide frequency:', file=f)
        for i in range(len(pepsList)):
            print(pepsList[i]+': '+str(pep_freq_tally[i])+' ('+str((int(10000*pep_freq_tally[i]/numseqs))/100)+'%)',file=f)
        print('\n'+str(num_freqs)+' most frequent alternate peptides found in sequences:',file=f)
        for i in range(len(pepsList)):
            if len(alts[i]) != 0:
                print(pepsList[i],file=f)
                these_alts = []
                for pep in alts[i]:
                    these_alts.append((pep,alts[i][pep]))
                print(str(sortAlts(these_alts,num_freqs))+'\n',file=f)
                #print(pepsList[i] + ' ' + str(alts[i])+'\n',file=f)
        print('\n*************************************************************\n',file=f)
        print(full_thing,file=f)
    return

#########################################################################################
#########################################################################################
#method to create spike sequences, again also takes parameter list from tkinter blocks
def getProtSeqs(fasta_path,paramsList):
    output_path = 'outputs/'+paramsList[3]
    search_range = paramsList[1]
    match_length = paramsList[2]
    reference = paramsList[0][:match_length]
    reference2 = paramsList[0][-match_length:]
    


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
                        prots_list.append(thisone)
                        break
                    thisone += prot[j]
        return prots_list

    #needleman wunsch / levenshtein. try fogsaa?
    #https://www.nature.com/articles/srep01746
    def similarity(s1, s2):
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
        return dp[len(dp) - 1][len(dp[0]) - 1]

    #find the best protein that matches reference
    def findTarget(genome,reference):
        maxProt = (0, '')
        for i in range(3):
            these_prots = findProteins(genome, i)
            for prot in these_prots:
                if (len(prot)>=(len(paramsList[0])*0.8) and len(prot)<=(len(paramsList[0])*1.2)):
                    score = similarity(prot[0:len(reference)], reference) + similarity(prot[-len(reference2):],reference2)
                    if score >= maxProt[0]:
                        maxProt = (score, prot)
        return maxProt


    startTime = time.time()
    bad_ones = []
    numSeqs = 0
    full_thing = ''
    #write final output here
    with open(output_path, 'w') as f:
        global today
        print('genome_to_spike spike output for file ' + fasta_path + '. Generated on '+ str(today) + '\n',file=f)
        print('parameters:\nsearch range = ' + str(search_range) +'\nmatch length = ' + str(match_length) +'\nreference:\n'+str(paramsList[0])+'\n',file=f)

        with open(fasta_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                numSeqs += 1
                # find best protein
                bestProtein = findTarget(record.seq[search_range[0]:search_range[1]].upper(), reference)
                if len(bestProtein[1]) < (len(paramsList[0])*0.8) or len(bestProtein[1]) > (len(paramsList[0])*1.2) or bestProtein[0]<(len(reference)):
                    bad_ones.append((record.id,len(bestProtein)))
                else:
                    '''
                    print('>'+record.id, file=f)
                    print(bestProtein[1], file=f)
                    '''
                    full_thing += '\n>'+record.id+'\n'+bestProtein[1]
        #print all bad sequences at start of file
        print(str(len(bad_ones))+' sequence(s) where target protein could not be found (likely obscured start/stop codons):\n'+str(bad_ones),file=f)

        endTime = time.time()
        total = endTime-startTime
        print('\nExecution time: '+(str(total))+' seconds for '+str(numSeqs)+' genomes\n\n*********************************************************',file=f)
        print('\n'+full_thing,file=f)
    return 'hi'
#########################################################################################
#########################################################################################
def prot_to_report(fasta_path,protParamz):
    output_path = 'outputs/'+protParamz[3]
    enz = protParamz[2]
    pepsList = protParamz[0]
    dist_algo = protParamz[1]
    num_freqs = 5

    min_pep_length, max_pep_length = 1000, 0
    for i in range(len(pepsList)):
        if len(pepsList[i]) >= max_pep_length:
            max_pep_length = len(pepsList[i])
        if len(pepsList[i]) <= min_pep_length:
            min_pep_length = len(pepsList[i])

    min_pep_length -= 3
    max_pep_length += 3

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

    def levenshtein(s1, s2):
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
                    if alg == 'hamming':
                        this_dist = hamming(interest[i], pep)
                    if alg == 'levenshtein':
                        this_dist = levenshtein(interest[i], pep)
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
            if enz == 'trypsin':
                pepper = tryptic(record.seq)
            if enz == 'alp':
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

    endTime = time.time()
    total_time = endTime - startTime

    with open(output_path, 'w') as f:
        print('peptide report for file ' + fasta_path + ' generated on ' + today, file=f)
        print('\nExecution time: ' + (str(total_time)) + ' seconds for ' + str(num_seqs) + ' sequences\n', file=f)
        for i in range(len(pepsList)):
            print(pepsList[i] + ': ' + str(pep_freq_tally[i]) + ', ' + str(100 * pep_freq_tally[i] / num_seqs)[
                                                                       0:5] + '%\n', file=f)

        print('\n' + str(num_freqs) + ' most frequent alternate peptides found in sequences:', file=f)
        for i in range(len(pepsList)):
            if len(alts[i]) != 0:
                print(pepsList[i], file=f)
                these_alts = []
                for pep in alts[i]:
                    these_alts.append((pep, alts[i][pep]))
                print(str(sortAlts(these_alts, num_freqs)) + '\n', file=f)

        print('\n*************************************************\n', file=f)
        print(details, file=f)
    return 





################################
    #view/edit profiles
################################
def viewEditProf():
    #read peptide profiles from saved text file to peptide_profile_dict
    with open('peptide_profiles.txt','r') as file:
        peptide_profile_dict = ast.literal_eval(file.read())
    
    win = Toplevel()
    win.geometry('240x190')
    win.wm_title("View/Edit Profiles")
    
    pl_area = scrolledtext.ScrolledText(win, wrap=WORD, 
                                      width=14, height=4, pady=2,padx=2)
    pl_area.grid(column=0, row=0)

    def changePl(*args):
        pl_area.delete('1.0',END)
        pl_area.insert(END, peptide_profile_dict[prof.get()])
    def saveProf():
        peptide_profile_dict[prof.get()] = pl_area.get('1.0', END).replace('\n','').replace(' ','')
        with open('peptide_profiles.txt','w') as f:
            print(peptide_profile_dict,file=f)
    def deleteProf():
        del peptide_profile_dict[prof.get()]
        with open('peptide_profiles.txt','w') as f:
            print(peptide_profile_dict,file=f)
        pep_prof_options = list(peptide_profile_dict.keys())
        w = OptionMenu(win,prof,*pep_prof_options,command=changePl)
        w.grid(row=1,column=0)
    def newProf():
        z = Toplevel()
        z.geometry('220x120')
        z.wm_title("new profile")
        saveNewBtn = Button(z, text="Save", command=lambda: [createProfile(),z.destroy()])
        saveNewBtn.grid(row=11, column=0)
        z_text = Entry(z)
        z_text.grid(row=1,column=0)
        z_lbl = Label(z,text='name')
        z_lbl.grid(row=1,column=1)
        pl_area = scrolledtext.ScrolledText(z, wrap=WORD, 
                                      width=14, height=4, pady=2,padx=2)
        pl_area.grid(column=0, row=2)
        np_lbl = Label(z,text='enter peptides')
        np_lbl.grid(row=2,column=1)
        def createProfile():
            peptide_profile_dict[z_text.get()] = pl_area.get("1.0", END).replace('\n','').replace(' ','')
            with open('peptide_profiles.txt','w') as f:
                print(peptide_profile_dict,file=f)
            pep_prof_options = list(peptide_profile_dict.keys())
            w = OptionMenu(win,prof,*pep_prof_options,command=changePl)
            w.grid(row=1,column=0)
 
        
    pp_lbl = Label(win,text='select user defined\npeptide profile')
    pp_lbl.grid(row=1,column=1)
    pep_prof_options = list(peptide_profile_dict.keys())
    prof = StringVar(win)
    prof.set('none')
    w = OptionMenu(win,prof,*pep_prof_options,command=changePl)
    w.grid(row=1,column=0)
    save_profile_btn = Button(win, text = 'save' ,fg = 'black',
                       command=saveProf)
    save_profile_btn.grid(row=2,column=0)
    del_profile_btn = Button(win, text = 'delete' ,fg = 'black',
                       command=deleteProf)
    del_profile_btn.grid(row=3,column=0)
    new_profile_btn = Button(win, text = 'new' ,fg = 'black',
                       command=newProf)
    new_profile_btn.grid(row=4,column=0)




vp_lbl = Label(main_display, text = '\n')
vp_lbl.grid(row=9,column=0)
view_profiles_btn = Button(main_display, text = 'view protein/peptide profiles' ,fg = 'black',
                       command=viewEditProf)
view_profiles_btn.grid(row=10,column=0)



if not validateDirectories():
    hide_all_widgets(main_display)
    label_1 = Label(main_display, text = '       Necessary files/directories not detected       ',fg = 'red')
    label_1.grid(row=0,column=0)
    label_2 = Label(main_display, text = '''Required directories:\nsequences\n(stores nucleotide/amino acid sequences)\noutputs\n(stores report outputs)\n\nRequired files:\nprotein_profiles.txt\n(stores user-defined protein sequences)\npeptide_profiles.txt\n(stores user-defined peptide lists)''')
    label_2.grid(row=1,column=0)
    label_3 = Label(main_display, text = '\nCreate these files/directories now?')
    label_3.grid(row=2,column=0)
    build_btn = Button(main_display, text = 'Create' ,fg = 'red',
                       command=lambda:[build_btn.destroy(),label_2.destroy(),
                                       label_3.destroy(),label_1.destroy(),buildDirectories(),
                                       show_all_widgets(main_display)])
    build_btn.grid(row=3,column=0)


main_display.mainloop()
