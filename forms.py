# AUTHOR: Shiri Almog , shirialmog1@gmail.com
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import StringField, SubmitField, SelectField, BooleanField,FormField
from wtforms.validators import DataRequired, Length

class rsIDform(FlaskForm):
    rsID = StringField("Option 1: Enter rsID",id='rsID',validators=[DataRequired()])
    submitrsID = SubmitField('Submit')

class fetch_coor(FlaskForm):
    genome=SelectField("Genome Assembly", choices=[('hg16','Human (hg16)'),('hg18','Human (hg18)'),('hg19','Human (hg19)'),('hg38','Human (hg38)'),
                       ('papAnu2','Baboon (papAnu2)'),('galGal2', 'Chicken (galGal2)'),('galGal3', 'Chicken (galGal3)'),('galGal4', 'Chicken (galGal4)'),('galGal5', 'Chicken (galGal5)'),('panTro1','Chimp (panTro1)'),
                        ('bosTau2','Cow (bosTau2)'),('bosTau3','Cow (bosTau3)'),('bosTau4','Cow (bosTau4)'),('bosTau7','Cow (bosTau7)'),('bosTau8','Cow (bosTau8)'),
                        ('macFas5','Crab-eating macaque (macFas5)'),('canFam1','Dog (canFam1)'),('canFam3','Dog (canFam3)'),('fr3','Fugu (fr3)'),
                        ('nomLeu3',  'Gibbon (nomLeu3)'),('equCab1','Horse (equCab1)'),('equCab2','Horse (equCab2)'),('anoCar2',' Lizard (anoCar2)'),
                        ('calJac3','Marmoset (calJac3)'),('oryLat2','Medaka (oryLat2)'),('mm10','Mouse (mm10)'),('mm7','Mouse (mm7)'),('mm8','Mouse (mm8)'),('mm9','Mouse (mm9)'),
                                                   ('monDom4','Opossum (monDom4)'),('monDom5','Opossum (monDom5)'),('susScr2','Pig (susScr2)'),('susScr3','Pig (susScr3)'),
                                                   ('ornAna1','Platypus (ornAna1)'),('ornAna2','Platypus (ornAna2)'),('oryCun2','Rabbit (oryCun2)'),('rn4','Rat (rn4)'),('rn5','Rat (rn5)'),('rn6','Rat (rn6)'),
                                                   ('rheMac2','Rhesus (rheMac2)'),('rheMac3','Rhesus (rheMac3)'),('rheMac8','Rhesus (rheMac8)'),('oviAri1','Sheep (oviAri1)'),('oviAri3','Sheep (oviAri3)'),
                                                   ('tetNig1','Tetraodon (tetNig1)'),('tetNig2','Tetraodon (tetNig2)'),('melGal1',' Turkey (melGal1)'),('taeGut1','Zebra finch (taeGut1)'),('taeGut2','Zebra finch (taeGut2)'),
                                                   ('danRer3','Zebrafish (danRer3)'),('danRer4','Zebrafish (danRer4)'),('danRer5','Zebrafish (danRer5)'),('danRer6','Zebrafish (danRer6)'),('danRer7','Zebrafish (danRer7)'),
                                                   ('danRer10','Zebrafish (danRer10)'),('sacCer1','S. cerevisiae (sacCer1)')])
    chromosome=StringField("Chromosome",validators=[DataRequired()])
    pos=StringField("Mutation position",validators=[DataRequired()])
    variation=SelectField('BE type',choices=[('A','A'),('C','C'),('G','G'),('T','T')],default='0',id="fromto")
    submitcoor = SubmitField('Submit')

class DNAform(FlaskForm):
    #rsID=StringField("If you know the rsID, enter it here:")
    PAM=StringField("PAM",id='PAM',validators=None)
    PAMstart = StringField("Activity window (distance from PAM)",id='start')
    PAMend = StringField("end",id='end')
    pambase=SelectField('BE type',choices=[('0','Select'),('1','CBE (C to T)'),('2','ABE (A to G)')],default='0',id="fromto")
    pamstream=SelectField('PAM orientation',choices=[('0','Select'),('U','Upstream'),('D','Downstream')],default='0',id="stream")

    upSeq = StringField("Upstream sequence:",id="Useq" ,validators=[DataRequired(), Length(min=25)])
    downSeq = StringField("Downstream sequence:", id="Dseq",
                          validators=[DataRequired(), Length(min=25)])
    mutation = StringField("Variant:",id="Mut", validators=[DataRequired(), Length(min=1, max=1)])
    WT = StringField("Reference:", id="Var",validators=[DataRequired(), Length(min=1, max=1)])
    readingFrame = StringField('Reading Frame:', id="RF", validators=None)
    #readingFrame = StringField('Reading Frame:',id="RF",validators=[DataRequired(),Length(min=1,max=1)])
    #showReadingFrame = BooleanField("Reading Frame:", default=True,validators=None,id='isRF')
    submit=SubmitField('Enter')
    #submitrsID = SubmitField('Submit')

class beffForm(FlaskForm):
    rsIDf=FormField(rsIDform)
    DNAf=FormField(DNAform)
    coorf=FormField(fetch_coor)

class contactForm(FlaskForm):
    message = StringField("Leave us a message here:", validators=[DataRequired()])
    submit = SubmitField('Enter')

class crisPAMform(FlaskForm):
    referencePAM = StringField("Enter reference Sequence here:", validators=[DataRequired(),Length(min=25)])
    variantPAM = StringField("Enter variant Sequence here:", validators=[DataRequired(),Length(min=25)])
    submit = SubmitField('Enter')

class uploadForm(FlaskForm):
    inputFile=FileField('Upload file here', validators=[DataRequired(), FileAllowed(['csv'])])
    submit=SubmitField('Upload')