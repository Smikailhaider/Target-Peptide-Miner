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
                if file == 'peptide_profiles.txt':
                    print("{'Test_profile':'peptide1,peptide2'}",file=f)
                if file == 'protein_profiles.txt':
                    print('>Test_protein\nMIKAIL',file=f)
                
                
    
   
#set up main window
main_display = Tk()
main_display.title('TPM')
main_display.geometry('255x255')

#create file selector
selected_file = ''
label_1 = Label(main_display, text = '         Select the sequences to be analyzed:         ')
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
        #genome_to_pep_btn['state'] = 'active'
        genome_to_prot_btn['state'] = 'active'
        prot_to_pep_btn['state'] = 'active'
button_explore = Button(main_display, text = "Browse Files",command = browseFiles)
button_explore.grid(row=1,column=0)

conf_lbl = Label(main_display, text = '[No file selected yet]\n')
conf_lbl.grid(column=0,row=2)

options_lbl = Label(main_display, text = 'Select desired operation:')
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
    final_lbl = Label(win, text='''An error occurred. Please ensure
you are using the correct sequences/files.''',fg='red')
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
    final_lbl.config(text='Successfully executed. File saved as:\n'+display_file)
    final_lbl.configure(font = ("Courier", 9, "bold"))
    final_lbl.grid(row=0,column=0)
    return
    

'''
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
    
    l = Label(win, text="Enter search parameters:")
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
    algo_options = ['hamming','NW']
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


    confirm_btn = Button(win, text="Confirm Parameters", command=lambda:[execReport()])
    confirm_btn.grid(row=8, column=0)

    def execReport():
        getPepParams = []
        getPepParams.append(pl_area.get('1.0', END).replace('\n', '').replace(' ','').split(','))
        getPepParams.append((int(sr_text.get()),int(sr2_text.get())))
        getPepParams.append(alg.get())
        getPepParams.append(report_area.get('1.0', END).replace('\n', '').replace(' ',''))
        try:
            genome_to_pep(selected_file,getPepParams)
            show_output_file(win,getPepParams[3])
        except:
            errorReset(win)
    return
    
genome_to_pep_btn = Button(main_display, text = 'genome -> peptide report' ,fg = 'red', command=genPepClick)
genome_to_pep_btn['state'] = 'disabled'
genome_to_pep_btn.grid(row=5,column=0)
'''




#driver for genome -> protein, popup window to specify method parameters
getProtParams = []
def getProt():
    
    #read protein profiles
    protein_profile_dict = {}
    protein_profile_dict['None'] = ''
    with open('protein_profiles.txt') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            protein_profile_dict[record.id] = record.seq
    
    
    win = Toplevel()
    win.geometry('270x290')
    win.wm_title("Input Parameters")
    l = Label(win, text="Enter search parameters:")
    l.grid(row=0, column=0)
    b = Button(win, text="Execute", command=lambda: [protExec()])
    b.grid(row=11, column=0)
    sr_text = Entry(win)
    sr_text.insert(0,'0')
    sr_text.grid(row=1,column=0)
    sr_start_lbl = Label(win,text='Search window start')
    sr_start_lbl.grid(row=1,column=1)

    sr2_text = Entry(win)
    sr2_text.insert(0,'0')
    sr2_text.grid(row=2,column=0)
    sr2_lbl = Label(win,text='Search window end')
    sr2_lbl.grid(row=2,column=1)

    tl_text = Entry(win)
    tl_text.insert(0,'75')
    tl_text.grid(row=3,column=0)
    tl_lbl = Label(win,text='Match length')
    tl_lbl.grid(row=3,column=1)

    rf_area = scrolledtext.ScrolledText(win, wrap=WORD, 
                                      width=14, height=4, pady=2,padx=2)
    rf_area.grid(column=0, row=4)
    rf_area.insert(END,'')
    rf_lbl = Label(win,text='Reference protein\nsequence')
    rf_lbl.grid(row=4,column=1)

    def changePl(*args):
        rf_area.delete('1.0',END)
        rf_area.insert(END, protein_profile_dict[prof.get()])
    
    prot_prof_user_selected = ''
    pp_lbl = Label(win,text='Select user defined\nprotein sequence')
    pp_lbl.grid(row=6,column=1)
    prot_prof_options = list(protein_profile_dict.keys())
    prof = StringVar(win)
    prof.set('None')
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
    report_lbl = Label(win,text='Output file \nname')
    report_lbl.grid(row=9,column=1)

    def protExec():
        getProtParams = []
        getProtParams.append(rf_area.get('1.0', END).replace('\n', '').replace(' ',''))
        getProtParams.append((int(sr_text.get()),int(sr2_text.get())))
        getProtParams.append(int(tl_text.get()))
        getProtParams.append(report_area.get('1.0', END).replace('\n', '').replace(' ',''))
        try:
            getProtSeqs(selected_file,getProtParams)
            show_output_file(win,getProtParams[3])
        except:
            errorReset(win)

genome_to_prot_btn = Button(main_display, text = 'Genome -> Target Protein' ,fg = 'red', command=getProt)
genome_to_prot_btn['state'] = 'disabled'
genome_to_prot_btn.grid(row=6,column=0)


#driver for protein -> peps, popup window to specify method parameters
getProtPepParams = []
def getProtReport():

    #***********************************************
    #read peptide profiles from saved text file to peptide_profile_dict
    with open('peptide_profiles.txt','r') as file:
        peptide_profile_dict = ast.literal_eval(file.read())
    peptide_profile_dict['None'] = ''
    #***********************************************
    
    win = Toplevel()
    win.geometry('270x310')
    win.wm_title("Input Parameters")
    l = Label(win, text="Enter report parameters")
    l.grid(row=0, column=0)

    tl_lbl = Label(win, text='Distance algorithm')
    tl_lbl.grid(row=3, column=1)
    algo_options = ['Hamming', 'NW']
    alg = StringVar(win)
    alg.set(algo_options[0])
    w = OptionMenu(win, alg, *algo_options)
    w.grid(row=3, column=0)

    tl_lbl2 = Label(win, text='Enzyme')
    tl_lbl2.grid(row=4, column=1)
    enz_options = ['Trypsin', 'A-LyticPrts']
    enz = StringVar(win)
    enz.set(enz_options[0])
    w2 = OptionMenu(win, enz, *enz_options)
    w2.grid(row=4, column=0)

    global pep_prof_user_selected
    pl_area = scrolledtext.ScrolledText(win, wrap=WORD,
                                        width=14, height=4, pady=2, padx=2)
    pl_area.grid(column=0, row=5)
    pl_lbl = Label(win, text='Enter peptides of\ninterest, separate with\ncommas')
    pl_lbl.grid(row=5, column=1)

    def changePl(*args):
        pl_area.delete('1.0',END)
        pl_area.insert(END, peptide_profile_dict[prof.get()])
    
    pep_prof_user_selected = ''
    pp_lbl = Label(win,text='Select user defined\npeptide profile')
    pp_lbl.grid(row=6,column=1)
    pep_prof_options = list(peptide_profile_dict.keys())
    prof = StringVar(win)
    prof.set('None')
    w2 = OptionMenu(win,prof,*pep_prof_options,command=changePl)
    w2.grid(row=6,column=0)
    ##########################
    pctExcl_text = Entry(win)
    pctExcl_text.insert(0,'-1')
    pctExcl_text.grid(row=7,column=0)
    pctExcl_lbl = Label(win, text='% threshold for\nalternate peptides')
    pctExcl_lbl.grid(row=7, column=1)

    b = Button(win, text="Execute",command=lambda: [prot_Exec()])
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
    report_lbl = Label(win,text='Report file\nname')
    report_lbl.grid(row=8,column=1)

    def prot_Exec():
        try:
            getProtPepParams = []
            getProtPepParams.append(pl_area.get('1.0', END).replace('\n', '').replace(' ','').split(','))
            getProtPepParams.append(alg.get())
            getProtPepParams.append(enz.get())
            getProtPepParams.append(report_area.get('1.0', END).replace('\n', '').replace(' ',''))
            getProtPepParams.append(pctExcl_text.get())
            prot_to_report(selected_file,getProtPepParams)
            show_output_file(win,getProtPepParams[3])
        
        except Exception as e:
            errorReset(win)
        
    return
prot_to_pep_btn = Button(main_display, text = 'Protein -> Peptide Report' ,fg = 'red', command=getProtReport)
prot_to_pep_btn['state'] = 'disabled'
prot_to_pep_btn.grid(row=7,column=0)


def limitChars(seq,limiter):
    output = ''
    for i in range(len(seq)):
        if i%limiter == 0 and i != 0:
            output += '\n'
        output += seq[i]
    return output


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
    #enter hamming or nw
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
                prot += codon_dict[genome[i:i + 3].replace('U','T')]
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
    '''
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
    '''
    def nw(s1, s2):
        rows, cols = (len(s1) + 1, len(s2) + 1)
        dp = [[(0,False) for i in range(cols)] for j in range(rows)]

        #affine gap penalties
        gap_open = -12
        gap_extend = -3
    
        for i in range(len(dp)):
            for j in range(len(dp[0])):
                #initialize array row/col 0
                if i == 0:
                    dp[i][j] = (-1*j,False)
                if j == 0:
                    dp[i][j] = (-1*i,False)
                
                if (i > 0 and j > 0):
                    '''
                    #get current residue match score
                    try:
                        match_score = pam_50_matrix[s1[i-1]][s2[j-1]]
                    except:
                        #handles unresolves bases/nucleotides: how forgiving should we be?
                        match_score = -1
                    '''
                    match_score = -3
                    if s1[i - 1] == s2[j - 1]:
                        match_score = 4
                    
                    match_res = dp[i-1][j-1][0] + match_score
                
                    #calculate (prospective) gap penalties
                    up_gap = dp[i-1][j][0]+gap_extend
                    if not dp[i-1][j][1]:
                        up_gap = dp[i-1][j][0]+gap_open
                    left_gap = dp[i][j-1][0]+gap_extend
                    if not dp[i-1][j][1]:
                        left_gap = dp[i][j-1][0]+gap_open
                    
                    #assign next value of 2d array
                    if (match_res >= up_gap) and (match_res >= left_gap):
                        dp[i][j] = (match_res,False)
                    elif (up_gap >= left_gap):
                        dp[i][j] = (up_gap,True)
                    else:
                        dp[i][j] = (left_gap,False)
        '''
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
        return len(s2)-dp[len(dp) - 1][len(dp[0]) - 1][0]
    

    #takes in a tryptic peptide list and one peptide of interest. returns best match with specified algorithm
    def bestMatch(trp,pep,full=False):
        if pep in trp:
            return 1
        min_dist = (len(pep) + 1,'')
        for match in trp:
            if match_algo == 'Hamming':
                thisDist = hamming(match,pep)
            if match_algo == 'NW':
                thisDist = nw(match,pep)
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
#method to create target protein sequences, again also takes parameter list from tkinter blocks
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

    #pam100,sc0.4
    pam_100_matrix = {'A': {'A': 3, 'R': -3, 'N': -1, 'D': -1, 'C': -3, 'Q': -1, 'E': 0, 'G': 0, 'H': -3, 'I': -1, 'L': -3, 'K': -3, 'M': -2, 'F': -4, 'P': 0, 'S': 1, 'T': 1, 'W': -6, 'Y': -4, 'V': 0}, 'R': {'A': -3, 'R': 6, 'N': -1, 'D': -3, 'C': -4, 'Q': 0, 'E': -3, 'G': -4, 'H': 1, 'I': -2, 'L': -4, 'K': 2, 'M': -1, 'F': -5, 'P': -1, 'S': -1, 'T': -2, 'W': 1, 'Y': -5, 'V': -3}, 'N': {'A': -1, 'R': -1, 'N': 4, 'D': 2, 'C': -5, 'Q': 0, 'E': 0, 'G': -1, 'H': 1, 'I': -2, 'L': -3, 'K': 1, 'M': -3, 'F': -4, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3}, 'D': {'A': -1, 'R': -3, 'N': 2, 'D': 5, 'C': -6, 'Q': 0, 'E': 3, 'G': -1, 'H': -1, 'I': -3, 'L': -5, 'K': -1, 'M': -4, 'F': -7, 'P': -3, 'S': -1, 'T': -1, 'W': -8, 'Y': -5, 'V': -3}, 'C': {'A': -3, 'R': -4, 'N': -5, 'D': -6, 'C': 8, 'Q': -7, 'E': -7, 'G': -4, 'H': -4, 'I': -3, 'L': -7, 'K': -7, 'M': -6, 'F': -6, 'P': -4, 'S': -1, 'T': -3, 'W': -8, 'Y': -1, 'V': -3}, 'Q': {'A': -1, 'R': 0, 'N': 0, 'D': 0, 'C': -7, 'Q': 5, 'E': 2, 'G': -3, 'H': 2, 'I': -3, 'L': -2, 'K': 0, 'M': -1, 'F': -6, 'P': 0, 'S': -2, 'T': -2, 'W': -6, 'Y': -5, 'V': -3}, 'E': {'A': 0, 'R': -3, 'N': 0, 'D': 3, 'C': -7, 'Q': 2, 'E': 5, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -7, 'P': -2, 'S': -1, 'T': -2, 'W': -8, 'Y': -5, 'V': -3}, 'G': {'A': 0, 'R': -4, 'N': -1, 'D': -1, 'C': -4, 'Q': -3, 'E': -1, 'G': 5, 'H': -4, 'I': -4, 'L': -5, 'K': -3, 'M': -4, 'F': -5, 'P': -2, 'S': 0, 'T': -2, 'W': -8, 'Y': -6, 'V': -2}, 'H': {'A': -3, 'R': 1, 'N': 1, 'D': -1, 'C': -4, 'Q': 2, 'E': -1, 'G': -4, 'H': 6, 'I': -4, 'L': -3, 'K': -2, 'M': -4, 'F': -3, 'P': -1, 'S': -2, 'T': -3, 'W': -3, 'Y': -1, 'V': -3}, 'I': {'A': -1, 'R': -2, 'N': -2, 'D': -3, 'C': -3, 'Q': -3, 'E': -3, 'G': -4, 'H': -4, 'I': 5, 'L': 1, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -3, 'T': 0, 'W': -6, 'Y': -2, 'V': 3}, 'L': {'A': -3, 'R': -4, 'N': -3, 'D': -5, 'C': -7, 'Q': -2, 'E': -4, 'G': -5, 'H': -3, 'I': 1, 'L': 5, 'K': -4, 'M': 2, 'F': 0, 'P': -3, 'S': -4, 'T': -3, 'W': -3, 'Y': -3, 'V': 0}, 'K': {'A': -3, 'R': 2, 'N': 1, 'D': -1, 'C': -7, 'Q': 0, 'E': -1, 'G': -3, 'H': -2, 'I': -3, 'L': -4, 'K': 5, 'M': 0, 'F': -6, 'P': -3, 'S': -1, 'T': -1, 'W': -5, 'Y': -5, 'V': -4}, 'M': {'A': -2, 'R': -1, 'N': -3, 'D': -4, 'C': -6, 'Q': -1, 'E': -3, 'G': -4, 'H': -4, 'I': 1, 'L': 2, 'K': 0, 'M': 7, 'F': -1, 'P': -3, 'S': -2, 'T': -1, 'W': -6, 'Y': -4, 'V': 1}, 'F': {'A': -4, 'R': -5, 'N': -4, 'D': -7, 'C': -6, 'Q': -6, 'E': -7, 'G': -5, 'H': -3, 'I': 0, 'L': 0, 'K': -6, 'M': -1, 'F': 7, 'P': -5, 'S': -3, 'T': -4, 'W': -1, 'Y': 4, 'V': -3}, 'P': {'A': 0, 'R': -1, 'N': -2, 'D': -3, 'C': -4, 'Q': 0, 'E': -2, 'G': -2, 'H': -1, 'I': -3, 'L': -3, 'K': -3, 'M': -3, 'F': -5, 'P': 6, 'S': 0, 'T': -1, 'W': -6, 'Y': -6, 'V': -2}, 'S': {'A': 1, 'R': -1, 'N': 1, 'D': -1, 'C': -1, 'Q': -2, 'E': -1, 'G': 0, 'H': -2, 'I': -3, 'L': -4, 'K': -1, 'M': -2, 'F': -3, 'P': 0, 'S': 3, 'T': 1, 'W': -2, 'Y': -3, 'V': -2}, 'T': {'A': 1, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': -2, 'H': -3, 'I': 0, 'L': -3, 'K': -1, 'M': -1, 'F': -4, 'P': -1, 'S': 1, 'T': 4, 'W': -6, 'Y': -3, 'V': 0}, 'W': {'A': -6, 'R': 1, 'N': -4, 'D': -8, 'C': -8, 'Q': -6, 'E': -8, 'G': -8, 'H': -3, 'I': -6, 'L': -3, 'K': -5, 'M': -6, 'F': -1, 'P': -6, 'S': -2, 'T': -6, 'W': 11, 'Y': -2, 'V': -7}, 'Y': {'A': -4, 'R': -5, 'N': -2, 'D': -5, 'C': -1, 'Q': -5, 'E': -5, 'G': -6, 'H': -1, 'I': -2, 'L': -3, 'K': -5, 'M': -4, 'F': 4, 'P': -6, 'S': -3, 'T': -3, 'W': -2, 'Y': 7, 'V': -3}, 'V': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -3, 'Q': -3, 'E': -3, 'G': -2, 'H': -3, 'I': 3, 'L': 0, 'K': -4, 'M': 1, 'F': -3, 'P': -2, 'S': -2, 'T': 0, 'W': -7, 'Y': -3, 'V': 5}}
    #pam50,sc0.5
    pam_50_matrix = {'A': {'A': 4, 'R': -4, 'N': -2, 'D': -1, 'C': -3, 'Q': -2, 'E': -1, 'G': 0, 'H': -4, 'I': -2, 'L': -3, 'K': -4, 'M': -2, 'F': -5, 'P': 0, 'S': 0, 'T': 0, 'W': -7, 'Y': -4, 'V': -1}, 'R': {'A': -4, 'R': 6, 'N': -3, 'D': -5, 'C': -4, 'Q': 0, 'E': -5, 'G': -5, 'H': 0, 'I': -3, 'L': -5, 'K': 1, 'M': -2, 'F': -5, 'P': -2, 'S': -1, 'T': -3, 'W': 0, 'Y': -6, 'V': -4}, 'N': {'A': -2, 'R': -3, 'N': 5, 'D': 2, 'C': -6, 'Q': -1, 'E': -1, 'G': -1, 'H': 1, 'I': -3, 'L': -4, 'K': 0, 'M': -4, 'F': -5, 'P': -3, 'S': 0, 'T': -1, 'W': -5, 'Y': -2, 'V': -4}, 'D': {'A': -1, 'R': -5, 'N': 2, 'D': 5, 'C': -8, 'Q': -1, 'E': 2, 'G': -1, 'H': -2, 'I': -4, 'L': -7, 'K': -2, 'M': -6, 'F': -8, 'P': -4, 'S': -2, 'T': -2, 'W': -8, 'Y': -6, 'V': -4}, 'C': {'A': -3, 'R': -4, 'N': -6, 'D': -8, 'C': 7, 'Q': -8, 'E': -8, 'G': -5, 'H': -4, 'I': -3, 'L': -8, 'K': -8, 'M': -7, 'F': -7, 'P': -4, 'S': -1, 'T': -4, 'W': -9, 'Y': -2, 'V': -3}, 'Q': {'A': -2, 'R': 0, 'N': -1, 'D': -1, 'C': -8, 'Q': 5, 'E': 1, 'G': -4, 'H': 1, 'I': -4, 'L': -3, 'K': -1, 'M': -2, 'F': -7, 'P': -1, 'S': -3, 'T': -3, 'W': -7, 'Y': -6, 'V': -4}, 'E': {'A': -1, 'R': -5, 'N': -1, 'D': 2, 'C': -8, 'Q': 1, 'E': 5, 'G': -2, 'H': -2, 'I': -3, 'L': -5, 'K': -2, 'M': -4, 'F': -8, 'P': -3, 'S': -2, 'T': -3, 'W': -9, 'Y': -5, 'V': -3}, 'G': {'A': 0, 'R': -5, 'N': -1, 'D': -1, 'C': -5, 'Q': -4, 'E': -2, 'G': 4, 'H': -5, 'I': -6, 'L': -6, 'K': -4, 'M': -5, 'F': -5, 'P': -3, 'S': 0, 'T': -3, 'W': -8, 'Y': -7, 'V': -3}, 'H': {'A': -4, 'R': 0, 'N': 1, 'D': -2, 'C': -4, 'Q': 1, 'E': -2, 'G': -5, 'H': 6, 'I': -5, 'L': -3, 'K': -3, 'M': -5, 'F': -3, 'P': -2, 'S': -3, 'T': -4, 'W': -4, 'Y': -1, 'V': -3}, 'I': {'A': -2, 'R': -3, 'N': -3, 'D': -4, 'C': -3, 'Q': -4, 'E': -3, 'G': -6, 'H': -5, 'I': 5, 'L': 0, 'K': -3, 'M': 0, 'F': -1, 'P': -5, 'S': -3, 'T': -1, 'W': -7, 'Y': -3, 'V': 2}, 'L': {'A': -3, 'R': -5, 'N': -4, 'D': -7, 'C': -8, 'Q': -3, 'E': -5, 'G': -6, 'H': -3, 'I': 0, 'L': 4, 'K': -4, 'M': 1, 'F': -1, 'P': -4, 'S': -5, 'T': -4, 'W': -3, 'Y': -4, 'V': -1}, 'K': {'A': -4, 'R': 1, 'N': 0, 'D': -2, 'C': -8, 'Q': -1, 'E': -2, 'G': -4, 'H': -3, 'I': -3, 'L': -4, 'K': 4, 'M': 0, 'F': -8, 'P': -3, 'S': -2, 'T': -1, 'W': -6, 'Y': -5, 'V': -5}, 'M': {'A': -2, 'R': -2, 'N': -4, 'D': -6, 'C': -7, 'Q': -2, 'E': -4, 'G': -5, 'H': -5, 'I': 0, 'L': 1, 'K': 0, 'M': 7, 'F': -2, 'P': -4, 'S': -3, 'T': -2, 'W': -7, 'Y': -6, 'V': 0}, 'F': {'A': -5, 'R': -5, 'N': -5, 'D': -8, 'C': -7, 'Q': -7, 'E': -8, 'G': -5, 'H': -3, 'I': -1, 'L': -1, 'K': -8, 'M': -2, 'F': 6, 'P': -6, 'S': -4, 'T': -5, 'W': -2, 'Y': 2, 'V': -4}, 'P': {'A': 0, 'R': -2, 'N': -3, 'D': -4, 'C': -4, 'Q': -1, 'E': -3, 'G': -3, 'H': -2, 'I': -5, 'L': -4, 'K': -3, 'M': -4, 'F': -6, 'P': 5, 'S': 0, 'T': -2, 'W': -8, 'Y': -7, 'V': -3}, 'S': {'A': 0, 'R': -1, 'N': 0, 'D': -2, 'C': -1, 'Q': -3, 'E': -2, 'G': 0, 'H': -3, 'I': -3, 'L': -5, 'K': -2, 'M': -3, 'F': -4, 'P': 0, 'S': 4, 'T': 1, 'W': -3, 'Y': -4, 'V': -3}, 'T': {'A': 0, 'R': -3, 'N': -1, 'D': -2, 'C': -4, 'Q': -3, 'E': -3, 'G': -3, 'H': -4, 'I': -1, 'L': -4, 'K': -1, 'M': -2, 'F': -5, 'P': -2, 'S': 1, 'T': 4, 'W': -7, 'Y': -3, 'V': -1}, 'W': {'A': -7, 'R': 0, 'N': -5, 'D': -8, 'C': -9, 'Q': -7, 'E': -9, 'G': -8, 'H': -4, 'I': -7, 'L': -3, 'K': -6, 'M': -7, 'F': -2, 'P': -8, 'S': -3, 'T': -7, 'W': 9, 'Y': -3, 'V': -9}, 'Y': {'A': -4, 'R': -6, 'N': -2, 'D': -6, 'C': -2, 'Q': -6, 'E': -5, 'G': -7, 'H': -1, 'I': -3, 'L': -4, 'K': -5, 'M': -6, 'F': 2, 'P': -7, 'S': -4, 'T': -3, 'W': -3, 'Y': 6, 'V': -4}, 'V': {'A': -1, 'R': -4, 'N': -4, 'D': -4, 'C': -3, 'Q': -4, 'E': -3, 'G': -3, 'H': -3, 'I': 2, 'L': -1, 'K': -5, 'M': 0, 'F': -4, 'P': -3, 'S': -3, 'T': -1, 'W': -9, 'Y': -4, 'V': 5}}
        
    
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

    #https://www.nature.com/articles/srep01746
    def similarity(s1, s2):
        rows, cols = (len(s1) + 1, len(s2) + 1)
        dp = [[(0,False) for i in range(cols)] for j in range(rows)]

        #affine gap penalties
        gap_open = -12
        gap_extend = -3
    
        for i in range(len(dp)):
            for j in range(len(dp[0])):
                #initialize array row/col 0
                if i == 0:
                    dp[i][j] = (-1*j,False)
                if j == 0:
                    dp[i][j] = (-1*i,False)
                
                if (i > 0 and j > 0):
                    
                    #get current residue match score
                    try:
                        match_score = pam_50_matrix[s1[i-1]][s2[j-1]]
                    except:
                        #handles unresolves bases/nucleotides: how forgiving should we be?
                        match_score = -3
                    '''
                    match_score = -3
                    if s1[i - 1] == s2[j - 1]:
                        match_score = 3
                    '''
                    match_res = dp[i-1][j-1][0] + match_score
                
                    #calculate (prospective) gap penalties
                    up_gap = dp[i-1][j][0]+gap_extend
                    if not dp[i-1][j][1]:
                        up_gap = dp[i-1][j][0]+gap_open
                    left_gap = dp[i][j-1][0]+gap_extend
                    if not dp[i-1][j][1]:
                        left_gap = dp[i][j-1][0]+gap_open
                    
                    #assign next value of 2d array
                    if (match_res >= up_gap) and (match_res >= left_gap):
                        dp[i][j] = (match_res,False)
                    elif (up_gap >= left_gap):
                        dp[i][j] = (up_gap,True)
                    else:
                        dp[i][j] = (left_gap,False)
        '''
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

    

    #find the best protein that matches reference
    def findTarget(genome,reference):
        maxProt = (-1, '')
        for i in range(3):
            these_prots = findProteins(genome, i)
            for prot in these_prots:
                if (len(prot)>=(len(paramsList[0])*0.97) and len(prot)<=(len(paramsList[0])*1.03)):
                    score_n = similarity(prot[0:len(reference)], reference)
                    score_c = similarity(prot[-len(reference2):],reference2)
                    if (score_n+score_c) >= maxProt[0]:
                        maxProt = (score_n+score_c, prot)
        return maxProt


    startTime = time.time()
    bad_ones = []
    numSeqs = 0
    full_thing = ''
    #write final output here
    with open(output_path, 'w') as f:
        global today
        print('Genome to protein output for file:\n' + fasta_path + '.\nGenerated on '+ str(today) + '.\n',file=f)
        print('Parameters:\nSearch range = ' + str(search_range) +'\nMatch length = ' + str(match_length) +'\nReference:\n'+limitChars(str(paramsList[0]),90)+'\n',file=f)

        with open(fasta_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                numSeqs += 1
                # find best protein
                if search_range == (0,0):
                    bestProtein = findTarget(record.seq.upper(), reference)
                else:
                    try:
                        bestProtein = findTarget(record.seq[search_range[0]:search_range[1]].upper(), reference)
                    except:
                        bestProtein = findTarget(record.seq[search_range[0]:len(record.seq)].upper(), reference)
                if len(bestProtein[1]) < (len(paramsList[0])*0.97) or len(bestProtein[1]) > (len(paramsList[0])*1.03) or bestProtein[0]<(len(reference)):
                    bad_ones.append((record.id,len(bestProtein)))
                else:
                    '''
                    print('>'+record.id, file=f)
                    print(bestProtein[1], file=f)
                    '''
                    full_thing += '\n>'+record.id+'\n'+limitChars(bestProtein[1],100)
        #print all bad sequences at start of file
        print(str(len(bad_ones))+' sequence(s) where target protein could not be found (potentially obscured start/stop codons):\n'+str(bad_ones),file=f)

        endTime = time.time()
        total = endTime-startTime
        print('\nExecution time: '+(str(round(total,2)))+' seconds for '+str(numSeqs)+' genomes\n\n*********************************************************',file=f)
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
    pctExclude = float(protParamz[4])

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
    total_time = endTime - startTime

    with open(output_path, 'w') as f:
        print('Peptide report for file: ' + fasta_path + ', generated on ' + today + '.', file=f)
        print('\nExecution time: ' + (str(round(total_time,2))) + ' seconds for ' + str(num_seqs) + ' sequences.\n', file=f)
        for i in range(len(pepsList)):
            print(pepsList[i] + ': ' + str(pep_freq_tally[i]) + ', ' + str(100 * pep_freq_tally[i] / num_seqs)[
                                                                       0:5] + '%\n', file=f)
        if pctExclude == -1:
            print('\n' + str(num_freqs) + ' most frequent alternate peptides found in sequences:', file=f)
        else:
            print('Most frequent alternate peptides found in sequences at >= '+(protParamz[4])+'% abundance:', file=f)
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
    return 





################################
    #view/edit peptide profiles
################################
def viewEditProf():
    #read peptide profiles from saved text file to peptide_profile_dict
    with open('peptide_profiles.txt','r') as file:
        try:
            peptide_profile_dict = ast.literal_eval(file.read())
        except:
            peptide_profile_dict = {}
    
    win = Toplevel()
    win.geometry('240x190')
    win.wm_title("View/Edit")
    
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
        if len(peptide_profile_dict) == 0:
            peptide_profile_dict['none'] = ''
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
view_profiles_btn = Button(main_display, text = 'View/Edit Peptide Profiles' ,fg = 'black',
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
