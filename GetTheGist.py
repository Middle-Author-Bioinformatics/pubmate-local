#!/usr/bin/env python3
import requests
import json
from datetime import datetime, timezone, timedelta
import re
from Bio import Entrez
from html import unescape
from openai import OpenAI
import argparse
import openai
import time
from collections import defaultdict
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.lib import colors
from PyPDF2 import PdfReader
from reportlab.platypus import Paragraph, SimpleDocTemplate, PageBreak, Spacer

apikey = open("/home/ark/MAB/bin/openai_apikey.txt", "r").read().strip()

client2 = OpenAI(api_key=apikey)

def save_to_pdf(results, output_file):
    """ Save PubMed results to a PDF file, ensuring each article starts on a new page """

    doc = SimpleDocTemplate(output_file, pagesize=letter, leftMargin=50, rightMargin=50, topMargin=50, bottomMargin=50)
    styles = getSampleStyleSheet()
    story = []

    for paper in results:
        # **Title of the Paper**
        story.append(Paragraph(f"<b>{paper['Title']}</b>", styles["Title"]))
        story.append(Spacer(1, 0.2 * inch))  # Add spacing

        # **Journal Name**
        story.append(Paragraph(f"<i>{paper['Journal']}</i>", styles["Normal"]))
        story.append(Spacer(1, 0.15 * inch))

        # **Publication Date**
        story.append(Paragraph(f"Publication Date: {paper['Publication Date']}", styles["Normal"]))
        story.append(Spacer(1, 0.15 * inch))

        # **PubMed & DOI Links**
        story.append(Paragraph(f"PubMed Link: {paper['PubMed Link']}", styles["Normal"]))
        story.append(Spacer(1, 0.15 * inch))
        story.append(Paragraph(f"DOI Link: {paper['DOI Link']}", styles["Normal"]))
        story.append(Spacer(1, 0.15 * inch))

        # **Authors**
        story.append(Paragraph(f"<b>Authors:</b> {paper['Authors']}", styles["Normal"]))
        story.append(Spacer(1, 0.2 * inch))

        # **Abstract**
        story.append(Paragraph(f"<b>Abstract:</b> {paper['Abstract']}", styles["Normal"]))
        story.append(Spacer(1, 0.2 * inch))

        # **Grant List & Keywords**
        story.append(Paragraph(f"<b>Funding:</b> {paper['Grant List']}", styles["Normal"]))
        story.append(Spacer(1, 0.15 * inch))
        story.append(Paragraph(f"<b>Keywords:</b> {paper['Keywords']}", styles["Normal"]))
        story.append(Spacer(1, 0.2 * inch))

        # **Page Break Between Articles**
        story.append(PageBreak())

    # **Build the PDF**
    doc.build(story)
    print(f"✅ Results saved to PDF file: {output_file}")


def format_gpt_text(text):
    """ Convert Markdown-style text to HTML for proper formatting in ReportLab """

    # Convert ### Headers to <b>Header</b>
    text = re.sub(r'(?m)^### (.*?)$', r'<b>\1</b>', text)

    # Convert **bold** to <b>bold</b>
    text = re.sub(r'\*\*(.*?)\*\*', r'<b>\1</b>', text)

    # Convert *italic* to <i>italic</i>
    text = re.sub(r'\*(.*?)\*', r'<i>\1</i>', text)

    return text

def chat_with_gpt(input_pdf, output_pdf, question, keywords):
    """ Calls GPT, reads input PDF (mandatory), and saves structured response to a separate output PDF """

    # **Step 1: Ensure Input PDF Exists and Extract Text**
    try:
        reader = PdfReader(input_pdf)
        pdf_text = "\n".join([page.extract_text() for page in reader.pages if page.extract_text()])
    except FileNotFoundError:
        print(f"❌ Error: Input PDF '{input_pdf}' not found. Cannot proceed.")
        return None  # Exit function if no input PDF is found
    except Exception as e:
        print(f"❌ Error reading input PDF: {e}")
        return None

    # **Step 2: Improved Prompt for Structured Output**
    structured_prompt = f"""
    You are analyzing a collection of research abstracts obtained from a PubMed search.
    The search was based on the following keywords: {keywords}.

    The extracted text from the PubMed abstracts is as follows:
    ---------------------------
    {pdf_text}
    ---------------------------

    Based on this document, please provide a structured response to the following question:
    "{question}"

    Format the response with:
    - A **summary of key findings** from the articles.
    - A **section on common themes or trends**.
    - A **section on gaps in knowledge or areas for future research**.
    - If applicable, include **bullet points or numbered lists** for clarity.
    - Ensure clear paragraph separation and professional formatting.
    - Please be clear which articles the findings are based on (i.e., cite the source).
    - Provide at least 1000 words

    Begin your response now:
    """

    # **Step 3: Call GPT to Get the Response**
    try:
        stream = client2.chat.completions.create(
            model="gpt-4o",
            messages=[{"role": "user", "content": structured_prompt}],
            stream=True
        )
        response_text = ""
        for chunk in stream:
            if chunk.choices[0].delta.content is not None:
                response_text += chunk.choices[0].delta.content

    except Exception as e:
        print(f"❌ GPT API Error: {e}")
        return None  # Exit function on GPT failure

    # **Step 4: Format the GPT Response**
    formatted_text = format_gpt_text(response_text)
    paragraphs = formatted_text.split("\n\n")  # Break text into paragraphs

    # **Step 5: Separate Bullet Points from Regular Text**
    bullet_points = []
    new_paragraphs = []

    for para in paragraphs:
        if para.startswith("- "):  # Identify bullet points
            bullet_points.append(para[2:].strip())  # Remove "- " and store separately
        else:
            if bullet_points:  # Add stored bullets before adding new paragraph
                new_paragraphs.append(bullet_points)
                bullet_points = []
            new_paragraphs.append(para.strip())

    if bullet_points:  # Append any remaining bullets at the end
        new_paragraphs.append(bullet_points)

    # **Step 6: Save GPT Response as a Structured PDF**
    doc = SimpleDocTemplate(output_pdf, pagesize=letter, leftMargin=50, rightMargin=50, topMargin=50, bottomMargin=50)
    styles = getSampleStyleSheet()
    story = []

    # **Add GPT response title**
    story.append(Paragraph("<b>GPT Response:</b>", styles["Title"]))
    story.append(Spacer(1, 0.2 * inch))

    # **Ensure each paragraph is formatted properly**
    for para in new_paragraphs:
        if isinstance(para, list):  # Handle bullet points separately
            for bullet in para:
                story.append(Paragraph(bullet, styles["BodyText"], bulletText="•"))  # Correct bullet usage
                story.append(Spacer(1, 0.1 * inch))
        else:
            story.append(Paragraph(para, styles["BodyText"]))
            story.append(Spacer(1, 0.2 * inch))

    # **Build the structured PDF**
    doc.build(story)
    print(f"✅ GPT response saved to structured PDF: {output_pdf}")

    return response_text  # Also return GPT's response for reference

def clean_text(text):
    # Remove HTML tags like <sub> and </sub>
    text = re.sub(r'<[^>]+>', '', text)

    # Unescape HTML entities like &kappa;, &ndash;, etc.
    text = unescape(text)

    return text


def search_pubmed(keywords, start_date, end_date, max_results, output_pdf_file1, output_pdf_file2, question):
    """ Fetch PubMed articles and save them all into one PDF file, including authors """
    Entrez.email = "your_email@example.com"  # Replace with your email

    query = ' AND '.join([f'({kw})' for kw in keywords])
    query = f'({query}) AND ("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication])'

    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="pub+date", retmode="xml")
    record = Entrez.read(handle)
    handle.close()

    id_list = record['IdList']
    pubmed_results = []

    if id_list:
        handle = Entrez.efetch(db="pubmed", id=",".join(id_list), rettype="medline", retmode="xml")
        papers = Entrez.read(handle)
        handle.close()

        for paper in papers['PubmedArticle']:
            title = clean_text(paper['MedlineCitation']['Article']['ArticleTitle'])
            pmid = paper['MedlineCitation']['PMID']
            journal = paper['MedlineCitation']['Article']['Journal']['Title']

            pub_date_info = paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
            pub_year = pub_date_info.get('Year', 'Unknown')
            pub_month = pub_date_info.get('Month', '')
            pub_day = pub_date_info.get('Day', '')
            publication_date = f"{pub_year} {pub_month} {pub_day}".strip()

            # ✅ Extract authors
            authors = []
            if 'AuthorList' in paper['MedlineCitation']['Article']:
                for author in paper['MedlineCitation']['Article']['AuthorList']:
                    last_name = author.get('LastName', '')
                    fore_name = author.get('ForeName', '')
                    authors.append(f"{fore_name} {last_name}" if fore_name else last_name)

            author_list = ", ".join(authors) if authors else "No authors listed"

            # ✅ Extract abstract (if available)
            abstract = None
            if 'Abstract' in paper['MedlineCitation']['Article']:
                abstract = " ".join(paper['MedlineCitation']['Article']['Abstract']['AbstractText'])

            # ✅ Extract DOI
            doi = None
            for identifier in paper['PubmedData']['ArticleIdList']:
                if identifier.attributes['IdType'] == 'doi':
                    doi = identifier

            # ✅ Extract grant information
            grant_list = []
            if 'GrantList' in paper['MedlineCitation']['Article']:
                for grant in paper['MedlineCitation']['Article']['GrantList']:
                    grant_id = grant.get('GrantID', 'N/A')
                    agency = grant.get('Agency', 'Unknown Agency')
                    grant_list.append(f"{agency} (Grant ID: {grant_id})")

            # ✅ Extract keywords
            keywords = []
            if 'KeywordList' in paper['MedlineCitation']:
                for keyword in paper['MedlineCitation']['KeywordList']:
                    keywords.extend(keyword)

            pubmed_results.append({
                "Title": title,
                "Journal": journal,
                "Publication Date": publication_date,
                "PubMed Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "DOI Link": f"https://doi.org/{doi}" if doi else "DOI not available",
                "Abstract": abstract if abstract else "No abstract available",
                "Grant List": ", ".join(grant_list) if grant_list else "No funding info",
                "Keywords": ", ".join(keywords) if keywords else "No keywords available",
                "Authors": author_list
            })

        # ✅ Save the results to the first PDF
        if pubmed_results:
            save_to_pdf(pubmed_results, output_pdf_file1)
            # print(f"✅ PDF saved successfully: {output_pdf_file1}")

            # ✅ Call chat_with_gpt() to generate a second analysis PDF
            chat_with_gpt(output_pdf_file1, output_pdf_file2, question, ",".join(keywords))
        else:
            print("⚠️ No search results found. No PDF generated.")
    else:
        print("No papers found for the specified criteria.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Upload files or directories to S3')
    parser.add_argument('--keywords', type=str, required=True, help='keywords to search, comma-separated')
    parser.add_argument('--question', type=str, required=False, help='question you are interested in (in quotes)')
    parser.add_argument('--output1', type=str, required=True, help='output first PDF file with results')
    parser.add_argument('--output2', type=str, required=True, help='output second PDF file with results')
    args = parser.parse_args()

    # Example usage
    keywords = args.keywords
    if args.question:
        question = args.question
    else:
        question = "What are the key findings of the articles"
    keywords = keywords.split(",")  # Keywords to search for in PubMed
    start_date = (datetime.now() - timedelta(days=100000)).strftime("%Y/%m/%d")
    end_date = datetime.now().strftime("%Y/%m/%d")
    max_results = 100  # Maximum number of results to fetch

    output_pdf_file1 = args.output1
    output_pdf_file2 = args.output2
    search_pubmed(keywords, start_date, end_date, max_results, output_pdf_file1, output_pdf_file2, question)
