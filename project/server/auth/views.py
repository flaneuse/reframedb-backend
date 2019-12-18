# project/server/auth/views.py


from flask import Blueprint, request, make_response, jsonify, session, url_for, render_template
from flask.views import MethodView

from data.example_data import example_data
from cdk_pywrapper.cdk_pywrapper import Compound

from project.server import app, bcrypt, db
from project.server.models import User, BlacklistToken
from project.server.token import generate_confirmation_token, confirm_token
from project.server.email import send_email

import pandas as pd
import wikidataintegrator as wdi
import requests
import os
import datetime

from elasticsearch import Elasticsearch
from elasticsearch.exceptions import RequestError
es = Elasticsearch()


auth_blueprint = Blueprint('auth', __name__)

data_dir = os.getenv('DATA_DIR')

assay_descrip = pd.read_csv(os.path.join(data_dir, '20180222_assay_descriptions.csv'), header=0)

plot_data = pd.read_csv(os.path.join(data_dir, '20180222_EC50_DATA_RFM_IDs_cpy.csv'), header=0)

#
# assay_data = pd.read_csv(data_dir + 'reframe_short_20170822.csv')
# gvk_dt = pd.read_csv(data_dir + 'gvk_data_to_release.csv')
# integrity_dt = pd.read_csv(data_dir + 'integrity_annot_20171220.csv')
#
# informa_dt = pd.read_csv(data_dir + 'informa_annot_20171220.csv')
#

# construct a map of InChI keys to Wikidata IDs
# ikey_wd_map = wdi.wdi_helpers.id_mapper('P235')
# wd_ikey_map = dict(zip(ikey_wd_map.values(), ikey_wd_map.keys()))
#
# print('wd ikey map length:', len(wd_ikey_map))


# retrieve chemical fingerprints and store them in a list, also store the associated IDs in a list
fingerprint_list = list()
id_list = list()


def calculate_tanimoto(fp_1, fp_2):
    intersct = fp_1.intersection(fp_2)
    return len(intersct)/(len(fp_1) + len(fp_2) - len(intersct))


def process_batch(batch):
    for count, hit in enumerate(batch['hits']['hits']):
        compound_id = hit['_id']
        id_list.append(compound_id)
        fingerprint = set(hit['_source']['fingerprint'])
        fingerprint_list.append(fingerprint)


body = {
    "_source": {
        "includes": ["fingerprint"],
    },
    "query": {
        "query_string": {

            "query": "*",
            "fields": ['fingerprint']
        }
    },
    "sort": [
        "_doc"
    ]
}

es.indices.refresh(index="reframe")

r = es.search(index="reframe", body=body, scroll='1m')
scroll_id = r['_scroll_id']
process_batch(r)

counter = 0
while True:
    r = es.scroll(scroll_id=scroll_id, scroll='1m')
    counter += 1
    #     print('batch processed:', counter)
    if len(r['hits']['hits']) == 0:
        break

    process_batch(r)

print('Total count of compound fingerprints:', len(fingerprint_list))

#
#
# def get_assay_data(qid):
#     tmp_dt = assay_data.loc[assay_data['wikidata'] == qid, :]
#
#     ad = list()
#
#     for c, x in tmp_dt.iterrows():
#         tmp_obj = dict()
#
#         # for k in x.keys():
#         #     tmp_obj.update({k: x[k]})
#
#         # only return the data really necessary for being rendered
#         datamode = x['datamode']
#
#         if datamode not in ['DECREASING', 'INCREASING', 'SUPER_ACTIVE']:
#             continue
#         elif datamode == 'DECREASING':
#             tmp_obj.update({'activity_type': 'IC50'})
#         elif datamode == 'INCREASING':
#             tmp_obj.update({'activity_type': 'EC50'})
#         elif datamode == 'SUPER_ACTIVE':
#             tmp_obj.update({'activity_type': 'SUPER ACTIVE'})
#
#         tmp_obj.update({'ac50': round(x['ac50'], 10)})
#         tmp_obj.update({'assay_title': x['assay_title']})
#         tmp_obj.update({'smiles': x['smiles']})
#         tmp_obj.update({'PubChem CID': str(x['PubChem CID'])})
#         tmp_obj.update({'wikidata': x['wikidata']})
#         tmp_obj.update({'calibr_id': x['calibr_id']})
#         tmp_obj.update({'inchi_key': x['inchi_key']})
#         tmp_obj.update({'ref': 'Calibr'})
#
#         ad.append(tmp_obj)
#
#     return ad
#
#
# def get_gvk_data(qid):
#     ikey = wd_ikey_map[qid]
#     print(ikey)
#
#     ad = list()
#
#     for c, x in gvk_dt.loc[gvk_dt['ikey'] == ikey, :].iterrows():
#         tmp_obj = {
#             'drug_name': x['drug_name'],
#             'phase': [{'label': y, 'qid': '', 'ref': 'GVK'} for y in x['phase'].split('; ')] if pd.notnull(x['phase']) else [],
#             'drug_roa': [{'label': y, 'qid': '', 'ref': 'GVK'} for y in x['drug_roa'].split('; ')] if pd.notnull(x['drug_roa']) else [],
#             'category': [{'label': y, 'qid': '', 'ref': 'GVK'} for y in x['category'].split('; ')] if pd.notnull(x['category']) else [],
#             'mechanism': [{'label': y, 'qid': '', 'ref': 'GVK'} for y in x['mechanism'].split('; ')] if pd.notnull(x['mechanism']) else [],
#             'synonyms': [{'label': y, 'qid': '', 'ref': 'GVK'} for y in x['synonyms'].split('; ')] if pd.notnull(x['synonyms']) else [],
#             'sub_smiles': x['sub_smiles'],
#         }
#
#         for cc, i in integrity_dt.loc[integrity_dt['ikey'] == ikey, :].iterrows():
#             if pd.notnull(i['status']):
#                 tmp_obj['phase'].extend(
#                     [{'label': y, 'qid': '', 'ref': 'Integrity'} for y in i['status'].split('; ')])
#             if pd.notnull(i['int_thera_group']):
#                 tmp_obj['category'].extend(
#                     [{'label': y, 'qid': '', 'ref': 'Integrity'} for y in i['int_thera_group'].split('; ')])
#             if pd.notnull(i['int_MoA']):
#                 tmp_obj['mechanism'].extend(
#                     [{'label': y, 'qid': '', 'ref': 'Integrity'} for y in i['int_MoA'].split('; ')])
#
#         for cc, i in informa_dt.loc[informa_dt['ikey'] == ikey, :].iterrows():
#             if pd.notnull(i['Global Status']):
#                 tmp_obj['phase'].extend(
#                     [{'label': y, 'qid': '', 'ref': 'Informa'} for y in i['Global Status'].split('; ')])
#             # if pd.notnull(i['int_thera_group']):
#                 # tmp_obj['category'].extend(
#                 #     [{'label': y, 'qid': '', 'ref': 'Informa'} for y in i['int_thera_group'].split('; ')])
#             if pd.notnull(i['Mechanism Of Action']):
#                 tmp_obj['mechanism'].extend(
#                     [{'label': y, 'qid': '', 'ref': 'Informa'} for y in i['Mechanism Of Action'].split('; ')])
#
#         ad.append(tmp_obj)
#
#     print(ad)
#
#     return ad


# --- LDH ---
def get_assay_details(assay_id):
    organisms = ['C. parvum', 'C. hominis', 'M. tuberculosis', 'Wolbachia', 'P. falciparum',
                 'Cryptosporidium', 'in vitro', 'in vivo', 'Mycobacterium tuberculosis', 'M. smegmatis']
    # def
    filtered = assay_descrip[assay_descrip.assay_id == assay_id].reset_index()

    # format the text
    filtered.assay_type = filtered.assay_type.apply(lambda x: x.lower())
    filtered.detection_method = filtered.detection_method.apply(lambda x: x.lower())
    filtered.drug_conc = filtered.drug_conc.apply(lambda x: x.replace('uM', ' µM'))

    return filtered


def get_assay_list():
    assays = []

    for idx, assay in assay_descrip.sort_values(['indication', 'assay_type', 'title_short']).iterrows():
        temp = {
            'assay_id': assay['assay_id'],
            'title_short': assay['title_short'],
            'indication': assay['indication'],
            'assay_type': assay['assay_type'],
            'summary':  assay['summary']
        }

        assays.append(temp)

    return assays


def get_dotplot_data(aid):
    def find_type(datamode):
        if datamode.lower() == 'decreasing':
            return 'IC'
        elif datamode.lower() == 'increasing':
            return 'EC'
        else:
            return 'unknown'

    def find_name(row):
        if pd.isnull(row.pubchem_label):
            return row.ID
        else:
            return row.pubchem_label

    assay_data = []
    # filter out data: selected assay, valid AC50 value
    filtered = plot_data.copy()[(plot_data['id'] == aid) & (plot_data['ac50'])]

    filtered['assay_type'] = filtered.datamode.apply(find_type)
    filtered = filtered.loc[filtered['assay_type'] != 'unknown'] # remove weird data modes

    filtered['url'] = filtered[['ikey', 'wikidata']].apply(
        lambda x: '/compound_data/{}'.format(x['ikey'])
        if pd.isnull(x['wikidata']) else '/compound_data/{};qid={}'.format(x['ikey'], x['wikidata']), axis=1)

    filtered['main_label'] = filtered.apply(find_name, axis=1)

    # group by unique ID and nest.
    # fncns = ['count', 'min', 'max', 'mean']
    # filtered = filtered.groupby(['ID']).agg(fncns)

    # sort values
    filtered = filtered.sort_values('ac50')

    # convert to the proper json-able structure
    for idx, cmpd in filtered.iterrows():
        temp = {
            'assay_title': cmpd.assay_title,
            'calibr_id': cmpd.calibr_id,
            'name': cmpd.main_label,
            'ac50': cmpd.ac50,
            'assay_type': cmpd.assay_type,
            'efficacy': cmpd.efficacy,
            'r_sq': cmpd.rsquared,
            'pubchem_id': cmpd['PubChem CID'] if pd.notnull(cmpd['PubChem CID']) else '',
            'url': cmpd.url
        }

        assay_data.append(temp)

    return assay_data


class RegisterAPI(MethodView):
    """
    User Registration Resource
    """

    def post(self):
        # get the post data
        post_data = request.get_json()

        # here, one needs to check with the Google ReCaptcha API whether ReCaptcha was sucessfully solved.
        # what would also be needed here is some kind of delay when a certain IP makes too many requests to either
        # signup oder login.

        recaptcha_token = post_data.get('recaptcha_token')
        params = {
            'secret': os.getenv('RECAPTCHA_SECRET_KEY'),
            'response': recaptcha_token
        }

        r = requests.post('https://www.google.com/recaptcha/api/siteverify', params=params).json()
        print(r)

        try:
            if not r['success']:
                response_object = {
                    'status': 'fail',
                    'message': 'ReCaptcha token could not be verified!'
                }
                return make_response(jsonify(response_object)), 401
        except KeyError:
            response_object = {
                'status': 'fail',
                'message': 'Some error occurred verifying ReCaptcha'
            }
            return make_response(jsonify(response_object)), 401

        # check if user already exists
        user = User.query.filter_by(email=post_data.get('email')).first()
        if not user:
            try:
                user = User(
                    email=post_data.get('email'),
                    password=post_data.get('password'),
                    confirmed=False
                )
                # insert the user
                db.session.add(user)
                db.session.commit()

                # generate email confirmation token
                token = generate_confirmation_token(user.email)
                confirm_url = app.config.get('FRONTEND_URL') + '#/confirm/' + token # url_for('auth.confirm_email', token=token, _external=True)
                html = render_template('activate.html', confirm_url=confirm_url)
                subject = "ReframeDB: Please confirm your email"
                send_email(user.email, subject, html)

                # generate the auth token
                #auth_token = user.encode_auth_token(user.id)
                responseObject = {
                    'status': 'success',
                    'message': 'A confirmation message has been sent to your email.',
                    #'auth_token': auth_token.decode()
                }
                return make_response(jsonify(responseObject)), 201
            except Exception as e:
                responseObject = {
                    'status': 'fail',
                    'message': 'Some error occurred. Please try again. ' + str(e)
                }
                return make_response(jsonify(responseObject)), 401
        else:
            responseObject = {
                'status': 'fail',
                'message': 'User already exists. Please Log in.',
            }
            return make_response(jsonify(responseObject)), 202


class LoginAPI(MethodView):
    """
    User Login Resource
    """
    def post(self):
        # get the post data
        post_data = request.get_json()
        try:
            # fetch the user data
            user = User.query.filter_by(
                email=post_data.get('email')
            ).first()
            if user:
                if bcrypt.check_password_hash(
                    user.password, post_data.get('password')
                ):
                    if user.confirmed:
                        auth_token = user.encode_auth_token(user.id)
                        if auth_token:
                            responseObject = {
                                'status': 'success',
                                'message': 'Successfully logged in.',
                                'auth_token': auth_token.decode(),
                                'confirmed': user.confirmed
                            }
                            return make_response(jsonify(responseObject)), 200
                    else:
                        responseObject = {
                            'status': 'fail',
                            'message': 'Please confirm your email before logging in.'
                        }
                        return make_response(jsonify(responseObject)), 500
                else:
                    responseObject = {
                        'status': 'fail',
                        'message': 'Incorrect password. Please try again.'
                    }
                    return make_response(jsonify(responseObject)), 500
            else:
                responseObject = {
                    'status': 'fail',
                    'message': 'User does not exist. Please Register.'
                }
                return make_response(jsonify(responseObject)), 404
        except Exception as e:
            print(e)
            responseObject = {
                'status': 'fail',
                'message': 'Try again'
            }
            return make_response(jsonify(responseObject)), 500


class UserAPI(MethodView):
    """
    User Resource
    """
    def get(self):
        # get the auth token
        auth_header = request.headers.get('Authorization')
        if auth_header:
            try:
                auth_token = auth_header.split(" ")[0]
            except IndexError:
                responseObject = {
                    'status': 'fail',
                    'message': 'Bearer token malformed.'
                }
                return make_response(jsonify(responseObject)), 401
        else:
            auth_token = ''
        if auth_token:
            resp = User.decode_auth_token(auth_token)
            if not isinstance(resp, str):
                user = User.query.filter_by(id=resp).first()
                responseObject = {
                    'status': 'success',
                    'data': {
                        'user_id': user.id,
                        'email': user.email,
                        'admin': user.admin,
                        'registered_on': user.registered_on
                    }
                }
                return make_response(jsonify(responseObject)), 200
            responseObject = {
                'status': 'fail',
                'message': resp
            }
            return make_response(jsonify(responseObject)), 401
        else:
            responseObject = {
                'status': 'fail',
                'message': 'Provide a valid auth token.'
            }
            return make_response(jsonify(responseObject)), 401


class LogoutAPI(MethodView):
    """
    Logout Resource
    """
    def post(self):
        # get auth token
        auth_header = request.headers.get('Authorization')
        if auth_header:
            auth_token = auth_header.split(" ")[0]
        else:
            auth_token = ''
        if auth_token:
            resp = User.decode_auth_token(auth_token)
            if not isinstance(resp, str):
                # mark the token as blacklisted
                blacklist_token = BlacklistToken(token=auth_token)
                try:
                    # insert the token
                    db.session.add(blacklist_token)
                    db.session.commit()
                    responseObject = {
                        'status': 'success',
                        'message': 'Successfully logged out.'
                    }
                    return make_response(jsonify(responseObject)), 200
                except Exception as e:
                    responseObject = {
                        'status': 'fail',
                        'message': e
                    }
                    return make_response(jsonify(responseObject)), 200
            else:
                responseObject = {
                    'status': 'fail',
                    'message': resp
                }
                return make_response(jsonify(responseObject)), 401
        else:
            responseObject = {
                'status': 'fail',
                'message': 'Provide a valid auth token.'
            }
            return make_response(jsonify(responseObject)), 403

class confirmEmail(MethodView):
    """
    Confirm Email Resource
    """
    def post(self):
        post_data = request.get_json()
        try:
            email = confirm_token(post_data.get('token'))
            if not email:
                responseObject = {
                    'status': 'fail',
                    'message': 'The confirmation link is invalid'
                }
                return make_response(jsonify(responseObject)), 500
        except:
            responseObject = {
                'status': 'fail',
                'message': 'The confirmation link is invalid or has expired'
            }
            return make_response(jsonify(responseObject)), 500
        user = User.query.filter_by(email=email).first_or_404()
        if user.confirmed:
            responseObject = {
                'status': 'fail',
                'message': 'Account already confirmed. Please login.'
            }
            return make_response(jsonify(responseObject)), 501
        else:
            user.confirmed = True
            user.confirmed_on = datetime.datetime.now()
            db.session.add(user)
            db.session.commit()
            responseObject = {
                'status': 'success',
                'message': 'You have confirmed your email. Please login.'
            }
            return make_response(jsonify(responseObject)), 200

class confirmEmailLink(MethodView):
    """
    Create new confirm link
    """
    def post(self):
        post_data = request.get_json()
        try:
            user = User.query.filter_by(
                email=post_data.get('email')
            ).first()
            if user:
                if not user.confirmed:
                    # generate email confirmation token
                    token = generate_confirmation_token(user.email)
                    confirm_url = app.config.get('FRONTEND_URL') + '#/confirm/' + token # url_for('auth.confirm_email', token=token, _external=True)
                    html = render_template('activate.html', confirm_url=confirm_url)
                    subject = "ReframeDB: Please confirm your email"
                    send_email(user.email, subject, html)

                    # generate the auth token
                    #auth_token = user.encode_auth_token(user.id)
                    responseObject = {
                        'status': 'success',
                        'message': 'A confirmation message has been sent to your email.',
                        #'auth_token': auth_token.decode()
                    }
                    return make_response(jsonify(responseObject)), 201
                else:
                    responseObject = {
                        'status': 'fail',
                        'message': 'You have already confirmed your email. Please login.'
                    }
                    return make_response(jsonify(responseObject)), 501
            else:
                responseObject = {
                    'status': 'fail',
                    'message': 'That user does not exist. Please register.'
                }
                return make_response(jsonify(responseObject)), 500
        except Exception as e:
            print(e)
            responseObject = {
                'status': 'fail',
                'message': 'Could not find user, please register or try again'
            }
            return make_response(jsonify(responseObject)), 500

class resetPassword(MethodView):
    """
    Reset the users password
    """
    def post(self):
        post_data = request.get_json()
        try:
            user = user = User.query.filter_by(id=post_data.get('user_id')).first()
            if user:
                user.password = bcrypt.generate_password_hash(
                    post_data.get('password'), app.config.get('BCRYPT_LOG_ROUNDS')
                ).decode()
                db.session.add(user)
                db.session.commit()
                responseObject = {
                    'status': 'success',
                    'message': 'You have changed your password. Please Login.'
                }
                return make_response(jsonify(responseObject)), 200
            else:
                responseObject = {
                    'status': 'fail',
                    'message': 'User could not be found.'
                }
                return make_response(jsonify(responseObject)), 501
        except Exception as e:
            print(e)
            responseObject = {
                'status': 'fail',
                'message': 'That user could not be found. Please register.'
            }
            return make_response(jsonify(responseObject)), 500

class resetPasswordCheck(MethodView):
    """
    Check if the reset link is valid
    """
    def post(self):
        post_data = request.get_json()
        try:
            email = confirm_token(post_data.get('token'))
            if not email:
                responseObject = {
                    'status': 'fail',
                    'message': 'The reset password link is invalid'
                }
                return make_response(jsonify(responseObject)), 500
        except:
            responseObject = {
                'status': 'fail',
                'message': 'The reset password link is invalid or has expired'
            }
            return make_response(jsonify(responseObject)), 500
        user = User.query.filter_by(email=email).first_or_404()
        if user.confirmed:
            responseObject = {
                'status': 'success',
                'message': 'Please enter your new password.',
                'data': {
                    'user_id': user.id,
                    'email': user.email
                }
            }
            return make_response(jsonify(responseObject)), 200
        else:
            # This response should never happen
            responseObject = {
                'status': 'fail',
                'message': 'You have not confirmed your email. Please confirm your email or create an account.'
            }
            return make_response(jsonify(responseObject)), 501

class resetPasswordLink(MethodView):
    """
    Create new reset password link
    """
    def post(self):
        post_data = request.get_json();
        try:
            user = User.query.filter_by(
                email=post_data.get('email')
            ).first()
            if user:
                if not user.confirmed:
                    responseObject = {
                        'status': 'fail',
                        'message': 'You have not confirmed your email. You cannot reset your password without a valid email.',
                        #'auth_token': auth_token.decode()
                    }
                    return make_response(jsonify(responseObject)), 500
                else:
                    # generate email confirmation token
                    token = generate_confirmation_token(user.email)
                    reset_pass_url = app.config.get('FRONTEND_URL') + '#/reset_pass/' + token # url_for('auth.confirm_email', token=token, _external=True)
                    html = render_template('password.html', reset_pass_url=reset_pass_url)
                    subject = "ReframeDB: Reset your password"
                    send_email(user.email, subject, html)

                    responseObject = {
                        'status': 'success',
                        'message': 'A new password link has been sent to your email.'
                    }
                    return make_response(jsonify(responseObject)), 200
            else:
                responseObject = {
                    'status': 'fail',
                    'message': 'That user does not exist. Please register.'
                }
                return make_response(jsonify(responseObject)), 500
        except Exception as e:
            print(e)
            responseObject = {
                'status': 'fail',
                'message': 'Could not find user, please register or try again'
            }
            return make_response(jsonify(responseObject)), 500

# class AssayDataAPI(MethodView):
#     """
#     Assaydata resource
#     """
#
#     def __init__(self):
#         pass
#
#     def get(self):
#         # get the auth token
#         auth_header = request.headers.get('Authorization')
#
#         print(auth_header)
#         args = request.args
#         qid = args['qid']
#
#         if auth_header:
#             try:
#                 auth_token = auth_header.split(" ")[0]
#             except IndexError:
#                 responseObject = {
#                     'status': 'fail',
#                     'message': 'Bearer token malformed.'
#                 }
#                 return make_response(jsonify(responseObject)), 401
#         else:
#             auth_token = ''
#         if auth_token:
#             resp = User.decode_auth_token(auth_token)
#             if not isinstance(resp, str):
#                 user = User.query.filter_by(id=resp).first()
#                 if user.id:
#                     responseObject = get_assay_data(qid)
#
#                     return make_response(jsonify(responseObject)), 200
#             responseObject = {
#                 'status': 'fail',
#                 'message': resp
#             }
#             return make_response(jsonify(responseObject)), 401
#         else:
#             responseObject = {
#                 'status': 'fail',
#                 'message': 'Provide a valid auth token.'
#             }
#             return make_response(jsonify(responseObject)), 401
#
#
# class GVKDataAPI(MethodView):
#     """
#     GVKData resource
#     """
#
#     def __init__(self):
#         pass
#
#     def get(self):
#         # get the auth token
#         auth_header = request.headers.get('Authorization')
#
#         print(auth_header)
#         args = request.args
#         qid = args['qid']
#
#         if auth_header:
#             try:
#                 auth_token = auth_header.split(" ")[0]
#             except IndexError:
#                 responseObject = {
#                     'status': 'fail',
#                     'message': 'Bearer token malformed.'
#                 }
#                 return make_response(jsonify(responseObject)), 401
#         else:
#             auth_token = ''
#         if auth_token:
#             resp = User.decode_auth_token(auth_token)
#             if not isinstance(resp, str):
#                 user = User.query.filter_by(id=resp).first()
#                 if user.id:
#                     responseObject = get_gvk_data(qid)
#
#                     return make_response(jsonify(responseObject)), 200
#             responseObject = {
#                 'status': 'fail',
#                 'message': resp
#             }
#             return make_response(jsonify(responseObject)), 401
#         else:
#             responseObject = {
#                 'status': 'fail',
#                 'message': 'Provide a valid auth token.'
#             }
#             return make_response(jsonify(responseObject)), 401


class SearchAPI(MethodView):
    """
    Search with elasticsearch
    """

    def __init__(self):
        pass

    @staticmethod
    def retrieve_document(compound_id, qid=''):
        search_result = {
            'id': '',  # InChI key or other
            'qid': '',  # if available
            'main_label': '',  # WHO INN or other
            'aliases': [],
            'assay_types': [],  # list of available assay types
            'tanimoto': 0,
            'similar_compounds': [],
            'reframeid': False,
            'smiles': '',
            'pubchem': '',
            'properties': [
                {'name': 'screening collection', 'tooltip': 'physical compound available in screening collection',
                 'value': False},
                {'name': 'assay hits', 'tooltip': 'compound responded in 1+ assays', 'value': False},
                {'name': 'Wikidata', 'tooltip': 'annotation data available from Wikidata', 'value': False},
                {'name': 'GVK', 'tooltip': 'annotation data available from GVK-Excelra GOSTAR', 'value': False},
                {'name': 'Integrity', 'tooltip': 'annotation data available from Clarivate Integrity', 'value': False},
                {'name': 'Citeline', 'tooltip': 'annotation data available from Citeline Pharmaprojects', 'value': False},
            ]
        }

        wd_main_label = ''
        who_inn = ''

        if qid:
            wd_r = es.get(index='wikidata', doc_type='compound', id=qid)
            wd_data = wd_r['_source']
            if 'en' in wd_data['labels']:
                wd_main_label = wd_data['labels']['en']['value']

            try:
                if 'P2275' in wd_data['claims']:
                    who_inn = wd_data['claims']['P2275'][0]['mainsnak']['datavalue']['value']['text']
            except KeyError:
                pass

            if who_inn:
                search_result['main_label'] = who_inn
            elif wd_main_label:
                search_result['main_label'] = wd_main_label

        r = es.get(index='reframe', doc_type='compound', id=compound_id)

        search_result['id'] = r['_id']
        data = r['_source']

        if 'qid' in data:
            qid = data['qid']
            search_result['qid'] = qid
            if qid:
                search_result['properties'][2]['value'] = True

        if 'reframe_id' in data and len(data['reframe_id']) > 0:
            search_result['reframeid'] = True
            search_result['properties'][0]['value'] = True
        else:
            search_result['reframeid'] = False

        aliases = set()
        for vendor, i in [('gvk', 3), ('informa', 5), ('integrity', 4)]:
            if len(data[vendor]) == 0:
                continue

            search_result['properties'][i]['value'] = True

            if not who_inn or not wd_main_label:
                search_result['main_label'] = data[vendor]['drug_name'][0]
                aliases.update(data[vendor]['drug_name'][1:])
            else:
                aliases.update(data[vendor]['drug_name'])
                aliases.discard(who_inn)
                aliases.discard(wd_main_label)

            if 'synonyms' in data[vendor]:
                aliases.update(data[vendor]['synonyms'])

            if 'PubChem CID' in data[vendor]:
                search_result['pubchem'] = data[vendor]['PubChem CID']

            if 'smiles' in data[vendor] and not search_result['smiles']:
                search_result['smiles'] = data[vendor]['smiles']

        if not search_result['main_label'] and len(aliases) > 0:
            search_result['main_label'] = aliases.pop()
        elif not search_result['main_label']:
            # as a last resort, set the ID as label
            search_result['main_label'] = r['_id']

        aliases.discard(search_result['main_label'])
        aliases.discard(search_result['main_label'].upper())
        aliases.discard(search_result['main_label'].lower())
        aliases.discard(search_result['main_label'].title())

        search_result['aliases'] = list(aliases)

        unique_assays = set()
        for assay in data['assay']:
            unique_assays.add((assay['assay_title'], assay['assay_id'] if 'assay_id' in assay else ''))
            search_result['properties'][1]['value'] = True

        search_result['assay_types'] = list(unique_assays)

        # get similar compounds from separate index
        if es.exists(index='similarity', doc_type='compound', id=compound_id):
            similarity_results = es.get(index='similarity', doc_type='compound', id=compound_id)

            if 'similar_compounds' in data:
                search_result['similar_compounds'] = similarity_results['_source']['similar_compounds']

        return search_result

    def get_compound_container(self, chem_id, suppress_hydrogens=True):
        try:
            return Compound(compound_string=chem_id, identifier_type='smiles', suppress_hydrogens=suppress_hydrogens)
        except ValueError as e:
            try:
                return Compound(compound_string=chem_id, identifier_type='inchi', suppress_hydrogens=suppress_hydrogens)
            except ValueError as e:
                raise e

    @staticmethod
    def exec_freetext_search(search_term, indices=('reframe', 'wikidata')):
        body = {
            "from": 0, "size": 100,
            "query": {
                "query_string": {
                    "default_operator": "AND",

                    "query": "{}*".format(search_term)
                }
            }
        }

        # body = {
        #     "from": 0, "size": 100,
        #     "query": {
        #         "multi_match": {
        #                 "query": "{}".format(search_term),
        #                 "operator": "AND",
        #
        #
        #
        #
        #         }
        #     }
        # }

        try:
            res = es.search(index=indices, body=body)
        except Exception as e:
            raise e

        results = []
        for x in res['hits']['hits']:
            compound_id = ''
            qid = ''

            if x['_index'] == 'reframe':
                compound_id = x['_id']
            elif x['_index'] == 'wikidata':
                qid = x['_id']
                # get InChI key
                compound_id = x['_source']['claims']['P235'][0]['mainsnak']['datavalue']['value']

                # make sure that Wikidata item gets used (labels, etc)
                for c, y in enumerate(results):
                    if y['id'] == compound_id:
                        results[c] = SearchAPI.retrieve_document(compound_id=compound_id, qid=qid)

            found = False
            for y in results:
                if y['id'] == compound_id:
                    found = True

            if not found:
                results.append(SearchAPI.retrieve_document(compound_id=compound_id, qid=qid))

        return results

    def get(self):
        # get the auth token
        auth_header = request.headers.get('Authorization')

        print(auth_header)
        args = request.args
        # TODO: add input checks here
        search_term = args['query']
        search_type = args['type']

        search_mode = ''
        if 'mode' in args:
            search_mode = args['mode']

        if 'tanimoto' in args and search_type == 'structure' and search_mode == 'similarity':
            tanimoto = float(args['tanimoto'])

            try:
                cmpnd = self.get_compound_container(search_term, suppress_hydrogens=True)
            except ValueError as e:
                return make_response(jsonify({
                    'status': 'fail',
                    'message': 'Tanimoto similarity cannot be calculated, no valid SMILES or InChI provided'
                })), 500

            searched_cmpnd_fp = {x for x in str(cmpnd.get_bitmap_fingerprint())[1:-1].split(', ')}

            found_cmpnds = []
            for c, x in enumerate(fingerprint_list):
                tnmt = calculate_tanimoto(searched_cmpnd_fp, x)
                if tnmt > tanimoto:
                    search_result = SearchAPI.retrieve_document(id_list[c])
                    search_result['tanimoto'] = tnmt
                    found_cmpnds.append(search_result)

            return make_response(jsonify({'status': 'success', 'results': found_cmpnds})), 200

        if search_type == 'string':

            analyze_body = {
                "analyzer": "standard",
                "text": search_term
            }

            analyze_result = es.indices.analyze(index='reframe', body=analyze_body)

            if len(analyze_result['tokens']) > 15:
                print('too complex search')
                return make_response(jsonify({'status': 'fail, too complex search query', 'results': []})), 400

            try:
                results = SearchAPI.exec_freetext_search(search_term)
            except Exception as e:
                return make_response(jsonify({'status': 'fail', 'results': []})), 500

            responseObject = {
                'status': 'success',
                'results': results
            }

            return make_response(jsonify(responseObject)), 200

        if search_type == 'structure' and (search_mode == 'exact' or search_mode == 'stereofree'):

            try:
                cmpnd = self.get_compound_container(search_term)
                ikey = cmpnd.get_inchi_key()
            except ValueError as e:
                return make_response(jsonify({
                    'status': 'fail',
                    'message': 'Similarity cannot be calculated, no valid SMILES or InChI provided'
                })), 500

            try:
                # cut off stere information from InChI key
                if search_mode == 'stereofree':
                    ikey = ikey[:15]

                results = SearchAPI.exec_freetext_search(ikey)
            except Exception as e:
                return make_response(jsonify({'status': 'fail', 'results': []})), 500

            responseObject = {
                'status': 'success',
                'results': results
            }

            return make_response(jsonify(responseObject)), 200

        return make_response(jsonify({
            'status': 'fail',
            'message': 'Malformed query'
        })), 500



        if auth_header:
            try:
                auth_token = auth_header.split(" ")[0]
            except IndexError:
                responseObject = {
                    'status': 'fail',
                    'message': 'Bearer token malformed.'
                }
                return make_response(jsonify(responseObject)), 401
        else:
            auth_token = ''
        if auth_token:
            resp = User.decode_auth_token(auth_token)
            if not isinstance(resp, str):
                user = User.query.filter_by(id=resp).first()
                if user.id:
                    responseObject = get_gvk_data(qid)

                    return make_response(jsonify(responseObject)), 200
            responseObject = {
                'status': 'fail',
                'message': resp
            }
            return make_response(jsonify(responseObject)), 401
        else:
            responseObject = {
                'status': 'fail',
                'message': 'Provide a valid auth token.'
            }
            return make_response(jsonify(responseObject)), 401


class MolfileAPI(MethodView):
    def __init__(self):
        pass

    def get(self):
        args = request.args
        compound_structure = args['compound_structure']

        try:
            compound = Compound(compound_string=compound_structure, identifier_type='smiles')
        except Exception as e:
            try:
                compound = Compound(compound_string=compound_structure, identifier_type='inchi')
                compound.get_inchi_key()
            except Exception as e:
                print(e)
                response = {
                    'status': 'fail',
                    'message': 'Invalid SMILES or InChI'
                }
                return make_response(jsonify(response)), 400

        molfile = compound.get_molfile()

        response = {
            'status': 'success',
            'molfile': molfile
        }

        return make_response(jsonify(response)), 200


class CompoundSVGAPI(MethodView):
    def __init__(self):
        pass

    def get(self):
        args = request.args
        compound_structure = args['structure']

        if 'format' in args and args['format'] == 'inchikey':
            search_results = SearchAPI.exec_freetext_search(compound_structure, indices=('reframe'))
            for x in search_results:
                if x['id'] == compound_structure:
                    compound_structure = x['smiles']
                    break

        try:
            compound = Compound(compound_string=compound_structure, identifier_type='smiles')
        except Exception as e:
            try:
                print('executed inchi conv')
                compound = Compound(compound_string=compound_structure, identifier_type='inchi')
                compound.get_inchi_key()
            except Exception as e:
                print(e)
                response = {
                    'status': 'fail',
                    'message': 'Invalid SMILES or InChI'
                }
                return make_response(jsonify(response)), 400

        svg_xml = compound.get_svg()

        response = {
            'status': 'success',
            'compound_svg': svg_xml
        }

        return make_response(jsonify(response)), 200


class DataAPI(MethodView):
    """
    Data resource
    """

    def __init__(self):
        pass

    def get(self):
        # get the auth token
        auth_header = request.headers.get('Authorization')

        print(auth_header)
        args = request.args
        qid = args['qid']
        # ikey = wd_ikey_map[qid]
        ikey = qid
        print(qid, ikey)

        if auth_header:
            try:
                auth_token = auth_header.split(" ")[0]
            except IndexError:
                responseObject = {
                    'status': 'fail',
                    'message': 'Bearer token malformed.'
                }
                return make_response(jsonify(responseObject)), 401
        else:
            auth_token = ''

        if auth_token:
            resp = User.decode_auth_token(auth_token)
            if not isinstance(resp, str):
                user = User.query.filter_by(id=resp).first()
                if user.id:
                    try:
                        res = es.get(index="reframe", doc_type='compound', id=ikey)
                    except Exception as e:
                        # if there is no Elasticsearch backend running, try to serve some example data.
                        if qid in example_data:
                            response_object = {qid: example_data[qid], "status": "success"}
                            return make_response(jsonify(response_object)), 200
                        else:
                            return make_response(jsonify({'status': 'fail', 'ikeys': []})), 500

                    responseObject = {
                        'status': 'success',
                        qid: res['_source']
                    }

                    return make_response(jsonify(responseObject)), 200

        else:
            responseObject = {
                'status': 'fail',
                'message': 'Provide a valid auth token.'
            }
            return make_response(jsonify(responseObject)), 401


# --- LDH ---
class AssayListAPI(MethodView):
    """
    Assay list resource
    """

    def __init__(self):
        pass

    def get(self):
        responseObject = get_assay_list()

        return make_response(jsonify(responseObject)), 200


# TODO: need a check if data returns something?
class AssayDetailsAPI(MethodView):
    """
    Assay description resource
    """

    def __init__(self):
        pass

    def get(self):
        args = request.args
        aid = args['aid']

        responseObject = get_assay_details(aid)
        return make_response(responseObject.to_json(orient="records")), 200


class PlotDataAPI(MethodView):
    """
    Assay list resource
    """

    def __init__(self):
        pass

    def get(self):
        # get the auth token
        auth_header = request.headers.get('Authorization')

        print(auth_header)
        args = request.args
        aid = args['aid']

        if auth_header:
            try:
                auth_token = auth_header.split(" ")[0]
            except IndexError:
                responseObject = {
                    'status': 'fail',
                    'message': 'Bearer token malformed.'
                }
                return make_response(jsonify(responseObject)), 401
        else:
            auth_token = ''
        if auth_token:
            resp = User.decode_auth_token(auth_token)
            if not isinstance(resp, str):
                user = User.query.filter_by(id=resp).first()
                if user.id:
                    responseObject = get_dotplot_data(aid)
                    print(responseObject)

# !! WARNING: jsonify can't handle NaN, since JSON hates NaN.

                    # responseObject = [{
                    #     "ac50": 1.25e-08,
                    #     "assay_title": "Crypto-C. parvum HCI proliferation assay - Wild Cp",
                    #     "assay_type": "IC",
                    #     "efficacy": -104.54743959999999,
                    #     "name": "drug891",
                    #     "pubchem_id": "CID2722",
                    #     "r_sq": 0.9641395209999999,
                    #     "test": "test",
                    #     "url": "/#/compound_data/Q239569"
                    # },
                    #     {
                    #     "ac50": 4.5199999999999994e-08,
                    #     "assay_title": "Crypto-C. parvum HCI proliferation assay - Wild Cp",
                    #     "assay_type": "IC",
                    #     "efficacy": -100.0,
                    #     "name": "drug888",
                    #     "pubchem_id": "CID2722",
                    #     "r_sq": 0.682526052,
                    #     "test": "test",
                    #     "url": "/#/compound_data/Q239569"
                    # }]

                    return make_response(jsonify(responseObject)), 200

            responseObject = {
                'status': 'fail',
                'message': resp
            }
            return make_response(jsonify(responseObject)), 401
        else:
            responseObject = {
                'status': 'fail',
                'message': 'Provide a valid auth token.'
            }
            return make_response(jsonify(responseObject)), 401


# define the API resources
registration_view = RegisterAPI.as_view('register_api')
login_view = LoginAPI.as_view('login_api')
user_view = UserAPI.as_view('user_api')
logout_view = LogoutAPI.as_view('logout_api')

# --- LDH ---
assay_list_view = AssayListAPI.as_view('assay_list_api')
assay_details_view = AssayDetailsAPI.as_view('assay_details_api')
plot_data_view = PlotDataAPI.as_view('plot_data_api')
# -----------


confirm_view = confirmEmail.as_view('confirm_email')
confirm_link = confirmEmailLink.as_view('confirm_link')
reset_pass_link = resetPasswordLink.as_view('reset_pass_link')
reset_pass_check = resetPasswordCheck.as_view('reset_pass_check')
reset_pass = resetPassword.as_view('reset_pass')
# assay_data_view = AssayDataAPI.as_view('assay_data_api')
# gvk_data_view = GVKDataAPI.as_view('gvk_data_api')
search_view = SearchAPI.as_view('search_api')
data_view = DataAPI.as_view('data_api')

molfile_view = MolfileAPI.as_view('molfile_conversion')
compound_svg_view = CompoundSVGAPI.as_view('compound_svg')

# add Rules for API Endpoints
auth_blueprint.add_url_rule(
    '/auth/register',
    view_func=registration_view,
    methods=['POST']
)
auth_blueprint.add_url_rule(
    '/auth/login',
    view_func=login_view,
    methods=['POST']
)
auth_blueprint.add_url_rule(
    '/auth/status',
    view_func=user_view,
    methods=['GET']
)
auth_blueprint.add_url_rule(
    '/auth/logout',
    view_func=logout_view,
    methods=['POST']
)
auth_blueprint.add_url_rule(
    '/auth/confirm',
    view_func=confirm_view,
    methods=['POST']
)
auth_blueprint.add_url_rule(
    '/auth/confirm/link',
    view_func=confirm_link,
    methods=['POST']
)
auth_blueprint.add_url_rule(
    '/auth/reset_pass',
    view_func=reset_pass,
    methods=['POST']
)
auth_blueprint.add_url_rule(
    '/auth/reset_pass/check',
    view_func=reset_pass_check,
    methods=['POST']
)
auth_blueprint.add_url_rule(
    '/auth/reset_pass/link',
    view_func=reset_pass_link,
    methods=['POST']
)
# auth_blueprint.add_url_rule(
#     '/assaydata',
#     view_func=assay_data_view,
#     methods=['GET'],
# )
# auth_blueprint.add_url_rule(
#     '/gvk_data',
#     view_func=gvk_data_view,
#     methods=['GET'],
# )
auth_blueprint.add_url_rule(
    '/search',
    view_func=search_view,
    methods=['GET'],
)
auth_blueprint.add_url_rule(
    '/data',
    view_func=data_view,
    methods=['GET'],
)
# --- assay lists (LDH) ---
auth_blueprint.add_url_rule(
    '/assay_list',
    view_func=assay_list_view,
    methods=['GET'],
)
auth_blueprint.add_url_rule(
    '/assay_details',
    view_func=assay_details_view,
    methods=['GET'],
)
auth_blueprint.add_url_rule(
    '/assaydata_plot',
    view_func=plot_data_view,
    methods=['GET'],
)

# chem API
auth_blueprint.add_url_rule(
    '/molfile',
    view_func=molfile_view,
    methods=['GET'],
)

auth_blueprint.add_url_rule(
    '/compound_svg',
    view_func=compound_svg_view,
    methods=['GET'],
)
