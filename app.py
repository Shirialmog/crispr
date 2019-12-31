from flask import Flask, render_template, url_for, flash, redirect,request,send_file, send_from_directory
from forms import DNAform, contactForm,crisPAMform, uploadForm,rsIDform,beffForm
from Scripts.BEsingle import MainBE
from Scripts.seqAnalysis import process_sequence
from Scripts.siteupload import Mainsiteupload
from Scripts.returncsv import returncsv
from Scripts.fetch_single_rsid import get_snp
from Scripts.fetch_genomic_seq import *
import os
import secrets

app = Flask(__name__)
app.config['SECRET_KEY']="3af987093f78a156771d3e00335bd60c"

@app.route('/', methods=['GET', 'POST'])
@app.route('/BE-FF', methods=['GET', 'POST'])
def input():
    form=beffForm()
    clean_dic=''
    quiet_dic=''
    refSeq=''
    mutSeq=''
    origMutSeq=''
    readingFrame=''
    locations_dic=''
    syn_quiet=''
    ## advanced ##
    pam=""
    pamstart=""
    pamend=""
    rsid_dic=''
    if "submit-rsid" in request.form and form.rsIDf.validate(form):
        rsid=form.rsIDf.rsID.data
        path=os.path.join(app.root_path, "user_files")
        tmp_dic=get_snp(rsid,path)
        rsid_dic=[]
        for num in range(len(tmp_dic)):
            tmp=[tmp_dic[num]['SNP'],tmp_dic[num]['CDS_flanking'],tmp_dic[num]['reading_frame']]
            if tmp not in rsid_dic:
                rsid_dic.append(tmp)


    if "submit-DNA" in request.form and form.DNAf.validate(form):
        readingFrame=form.DNAf.readingFrame.data
        if form.DNAf.PAM.data:
            pam=form.DNAf.PAM.data
        if form.DNAf.PAMstart.data:
            pamstart=int(form.DNAf.PAMstart.data)
        if form.DNAf.PAMend.data:
            pamend=int(form.DNAf.PAMend.data)
        #clean_dic, quiet_dic, refSeq, mutSeq, origMutSeq, locations_dic, syn_quiet = MainBE(form.upSeq.data,form.downSeq.data,
                                                                                            #form.mutation.data, form.WT.data, readingFrame)
        clean_dic,quiet_dic, refSeq,mutSeq, origMutSeq, locations_dic,syn_quiet = MainBE(form.DNAf.upSeq.data,form.DNAf.downSeq.data,form.DNAf.mutation.data,form.DNAf.WT.data,readingFrame,pam,pamstart,pamend,form.DNAf.pambase.data,form.DNAf.pamstream.data)
    return render_template('BE-FF.html', title="BE-FF", form=form, RF=readingFrame, result=(clean_dic,quiet_dic,refSeq,mutSeq,origMutSeq, locations_dic,syn_quiet,rsid_dic),rsdic=rsid_dic) # readingFrame=readingFrame, showReadingFrame=showReadingFrame)



@app.route('/help')
def help():
    return render_template('help.html', title='Help')
@app.route('/contact')
def contact():
    form=contactForm()
    return render_template('contact.html',title="About", form=form)

def save_file(input_field):
    random_hex=secrets.token_hex(8)
    f_name, f_ext=os.path.splitext(input_field.filename)
    filename=random_hex+f_ext
    filepath=os.path.join(app.root_path, "user_files",filename)
    input_field.save(filepath)
    return filepath,filename

@app.route('/sample')
def sample():
     filepath = os.path.join(app.root_path, "user_files", 'sample_file.csv')
     return send_file(filepath, attachment_filename='sample_file.csv', as_attachment=True)

@app.route('/download')
def download():
     filepath = os.path.join(app.root_path, "user_files", 'results.csv')
     return send_file(filepath, attachment_filename='results.csv', as_attachment=True)


@app.route('/upload',methods=['GET', 'POST'])
def upload():
    form=uploadForm()
    result=''
    filename=''
    showresults = False
    if form.validate_on_submit():
        filepath,filename=save_file(form.inputFile.data)
        matches, cleanMatchdic, quietMatchdic, rsltsDic=Mainsiteupload(filepath)
        returncsv(matches, cleanMatchdic, quietMatchdic,rsltsDic,filepath)
        showresults=True
        send_from_directory('user_files', filename)

    elif request.method=='GET':
        # showresults=False
        result='1'

    return render_template('upload.html',form=form,result=result,
                           filename=filename,showresults=showresults)

@app.route('/download2/<filename>', methods=['GET', 'POST'])
def download2(filename):
    return send_from_directory('user_files', filename)

if __name__ == '__main__':
    app.run(debug=True)
