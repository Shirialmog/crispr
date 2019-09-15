from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, SelectField, BooleanField
from wtforms.validators import DataRequired, Length

class DNAform(FlaskForm):
    upSeq = StringField("Upstream sequence:", validators=[DataRequired(), Length(min=25)])
    downSeq = StringField("Downstream sequence:",
                          validators=[DataRequired(), Length(min=25)])
    mutation = StringField("Variant:", validators=[DataRequired(), Length(min=1, max=1)])
    WT = StringField("Reference:", validators=[DataRequired(), Length(min=1, max=1)])
    readingFrame = SelectField("Reading Frame:", choices=[(1, 1), (2, 2), (3, 3)])
    showReadingFrame = BooleanField("Reading Frame:", default=False)
    submit=SubmitField('Enter')

class contactForm(FlaskForm):
    message = StringField("Leave us a message here:", validators=[DataRequired()])
    submit = SubmitField('Enter')

class crisPAMform(FlaskForm):
    referencePAM = StringField("Enter reference Sequence here:", validators=[DataRequired(),Length(min=25)])
    variantPAM = StringField("Enter variant Sequence here:", validators=[DataRequired(),Length(min=25)])
    submit = SubmitField('Enter')