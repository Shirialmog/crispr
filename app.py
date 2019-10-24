from flask import Flask, render_template, url_for, flash, redirect,request,send_file, send_from_directory
from forms import DNAform, contactForm,crisPAMform, uploadForm
from Scripts.BEsingle import MainBE
from Scripts.seqAnalysis import process_sequence
from Scripts.BEsiteupload import Mainupload
from Scripts.writeresults import write_results
import os
import secrets

app = Flask(__name__)
app.config['SECRET_KEY']="3af987093f78a156771d3e00335bd60c"
posts = [
    {
        'author': 'Shiri Almog',
        'title': 'Blog Post 1',
        'content': 'First post content',
        'date_posted': 'Sept 10, 2019'
    },
    {
        'author': '123',
        'title': 'Blog Post 2',
        'content': 'Second post content',
        'date_posted': 'April 21, 2018'
    }
]

@app.route('/')
@app.route('/home')
def hello_world():
    return render_template('home.html', posts=posts)
@app.route('/blog')
def hello_world2():
    return render_template('blog.html', posts=posts)

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/contact')
def contact():
    form=contactForm()
    return render_template('contact.html',title="contact", form=form)


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
     print (filepath)
     return send_file(filepath, attachment_filename='sample_file.csv', as_attachment=True)

@app.route('/download')
def download():
     filepath = os.path.join(app.root_path, "user_files", 'full_results.csv')
     return send_file(filepath, attachment_filename='full_results.csv', as_attachment=True)


@app.route('/upload',methods=['GET', 'POST'])
def upload():
    form=uploadForm()
    result=''
    filename=''
    showresults = False
    if form.validate_on_submit():
        filepath,filename=save_file(form.inputFile.data)
        matches, cleanMatchdic, quietMatchdic, rsltsDic=Mainupload(filepath)
        write_results(matches, cleanMatchdic, quietMatchdic,rsltsDic,'full_results.csv')
        showresults=True

    elif request.method=='GET':
        # showresults=False
        result='1'

    return render_template('upload.html',form=form,result=result,
                           filename='full_results.csv',showresults=showresults)

@app.route('/crisPAM', methods=['GET','POST'])
def crisPAM():
    form=crisPAMform()
    result=''
    if form.validate_on_submit():
        result=process_sequence(form.referencePAM.data,form.variantPAM.data)
    return render_template('crisPAM.html',title='crisPAM',form=form,result=result)


@app.route('/input', methods=['GET', 'POST'])
def input():
    form=DNAform()
    clean_dic=''
    quiet_dic=''
    refSeq=''
    mutSeq=''
    origMutSeq=''
    origRevSeq=''
    printRevCorSeq=''
    printPam=''
    readingFrame=''
    locations_dic=''
    if form.validate_on_submit():
        if form.showReadingFrame==False:
            readingFrame=0
        else:
            readingFrame=form.readingFrame.data
            print (readingFrame)
        clean_dic,quiet_dic, refSeq,mutSeq, origMutSeq, origRevSeq,locations_dic = MainBE(form.upSeq.data,form.downSeq.data,form.mutation.data,form.WT.data,readingFrame)
    return render_template('input.html', title="input", form=form, RF=readingFrame, result=(clean_dic,quiet_dic,refSeq,mutSeq,origMutSeq, origRevSeq,locations_dic)) # readingFrame=readingFrame, showReadingFrame=showReadingFrame)


@app.route('/download2/<filename>', methods=['GET', 'POST'])
def download2(filename):
    return send_from_directory('user_files', filename)

@app.route('/formTable', methods=['POST', 'GET'])
def formTable():
    if request.method=='POST':
        input1=request.form.getE
        input2 = request.form.table2
        input3 = request.form.table3
        input4 = request.form.table4
        return render_template('formTable.html', result=(input1,input2,input3,input4))

@app.route('/getTable', methods=['POST', 'GET'])
def getTable():
    result=0
    return render_template('getTable.html',result=result)

if __name__ == '__main__':
    app.run(debug=True)
