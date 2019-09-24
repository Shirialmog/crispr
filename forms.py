from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import StringField, SubmitField, SelectField, BooleanField
from wtforms.validators import DataRequired, Length

class DNAform(FlaskForm):
    rsID=StringField("If you know the rsID, enter it here:")
    upSeq = StringField("Upstream sequence:",id="Useq" ,validators=[DataRequired(), Length(min=25)])
    downSeq = StringField("Downstream sequence:", id="Dseq",
                          validators=[DataRequired(), Length(min=25)])
    mutation = StringField("Variant:",id="Var", validators=[DataRequired(), Length(min=1, max=1)])
    WT = StringField("Reference:", id="Mut",validators=[DataRequired(), Length(min=1, max=1)])
    readingFrame = StringField("Reading Frame:",id="RF",validators=None)
    showReadingFrame = BooleanField("Reading Frame:", default=True,validators=None)
    submit=SubmitField('Enter')
    submitrsID = SubmitField('Submit')

class contactForm(FlaskForm):
    message = StringField("Leave us a message here:", validators=[DataRequired()])
    submit = SubmitField('Enter')

class crisPAMform(FlaskForm):
    referencePAM = StringField("Enter reference Sequence here:", validators=[DataRequired(),Length(min=25)])
    variantPAM = StringField("Enter variant Sequence here:", validators=[DataRequired(),Length(min=25)])
    submit = SubmitField('Enter')

class uploadForm(FlaskForm):
    inputFile=FileField('Upload file here', validators=[FileAllowed(['csv'])])
    submit=SubmitField('Upload')