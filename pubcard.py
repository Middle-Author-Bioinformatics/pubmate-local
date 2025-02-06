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
import os
from collections import defaultdict
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.platypus import Image
from reportlab.platypus import Paragraph, SimpleDocTemplate, PageBreak, Spacer
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.lib import colors
from PyPDF2 import PdfReader
import base64

apikey = open("/home/ark/MAB/bin/openai_apikey.txt", "r").read().strip()
print(apikey)
client = OpenAI(api_key=apikey)

def generate_author_variations(last_name, first_name, middle_initial=""):
    """
    Generate different variations of the author's name based on the provided details.
    """
    name_variations = [
        rf"\b{last_name}\s+{first_name[0]}\.?{middle_initial}\b",  # Last FM
        rf"\b{first_name}\s+{middle_initial}\s*{last_name}\b",  # First M Last
        rf"\b{first_name}\s+{last_name}\b",  # First Last
        rf"\b{last_name}\s+{first_name[0]}\b"  # Last F
    ]
    if middle_initial:  # Add variations that include the middle initial optionally
        name_variations.append(rf"\b{last_name}\s+{first_name[0]}\.?\b")  # Last F (without middle initial)

    return name_variations


def generate_pubmed_author_query(last_name, first_name, middle_initial):
    """
    Generate multiple variations of author names for PubMed searches.
    If a middle initial is provided, it must be included in all variations.
    """
    first_initial = first_name[0]

    if middle_initial:  # Middle initial is required in all searches
        variations = [
            f'"{last_name} {first_initial}{middle_initial}"[Author]',  # Last FM
            f'"{first_name} {middle_initial} {last_name}"[Author]',  # First M Last
            f'"{last_name} {first_initial} {middle_initial}"[Author]',  # Last F M
        ]
    else:  # If no middle initial, allow variations without it
        variations = [
            f'"{last_name} {first_initial}"[Author]',  # Last F
            f'"{first_name} {last_name}"[Author]',  # First Last
            f'"{last_name} {first_initial}."[Author]',  # Last F. (with period)
        ]

    return " OR ".join(variations)  # Combine with OR for broader search


def generate_pubmed_author_query(last_name, first_name, middle_initial):
    """
    Generate multiple variations of author names for PubMed searches.
    If a middle initial is provided, it must be included in all variations.
    """
    first_initial = first_name[0]

    if middle_initial:  # Middle initial is required in all searches
        variations = [
            f'"{last_name} {first_initial}{middle_initial}"[Author]',  # Last FM
            f'"{first_name} {middle_initial} {last_name}"[Author]',  # First M Last
            f'"{last_name} {first_initial} {middle_initial}"[Author]',  # Last F M
        ]
    else:  # If no middle initial, allow variations without it
        variations = [
            f'"{last_name} {first_initial}"[Author]',  # Last F
            f'"{first_name} {last_name}"[Author]',  # First Last
            f'"{last_name} {first_initial}."[Author]',  # Last F. (with period)
        ]

    return " OR ".join(variations)  # Combine with OR for broader search


def extract_abstracts_from_pdf(pdf_path, last_name, first_name, middle_initial=""):
    """
    Extracts abstracts from a given PDF file, filtering for cases where the author appears
    as the first, second, or last author.
    """
    reader = PdfReader(pdf_path)
    pdf_text = "\n".join([page.extract_text() for page in reader.pages if page.extract_text()])

    # Generate name variations dynamically
    author_patterns = generate_author_variations(last_name, first_name, middle_initial)

    # Extract abstracts and filter based on author position
    abstracts = []
    abstracts_found = re.findall(r"Abstract:\s*(.*?)\nFunding:", pdf_text, re.DOTALL)

    for abstract in abstracts_found:
        author_match = False
        # Extract author list
        author_list_match = re.search(r"Authors:\s*\[([^\]]+)\]", pdf_text)
        if author_list_match:
            author_list = author_list_match.group(1).split(", ")

            # Check if author is first, second, or last
            if len(author_list) > 1:
                first_author, second_author = author_list[0], author_list[1]
                last_author = author_list[-1]
                author_positions = [first_author, second_author, last_author]
            else:
                author_positions = [author_list[0]]

            # Check for match with any author name variation
            for author_pattern in author_patterns:
                if any(re.search(author_pattern, name, re.IGNORECASE) for name in author_positions):
                    author_match = True
                    break

        if author_match:
            abstracts.append(abstract.strip())

    return abstracts

def summarize_in_batches(abstracts, batch_size=3):
    """ Summarizes abstracts in smaller batches to avoid token limits. """
    summary = ''
    abstract_list = abstracts.split("\n\n")  # Split abstracts into individual items

    for i in range(0, len(abstract_list), batch_size):
        batch = "\n\n".join(abstract_list[i:i+batch_size])  # Take a few abstracts at a time

        try:
            response = client.chat.completions.create(
                model="gpt-4",
                messages=[
                    {"role": "system", "content": "Summarize the following abstracts into key research themes."},
                    {"role": "user", "content": batch}
                ],
                max_tokens=5000  # Keep responses short
            )
            summaries.append(response.choices[0].message.content.strip())

        except Exception as e:
            print(f"⚠️ Error summarizing batch: {e}")

    return "\n".join(summaries) if summaries else None

def generate_document_image(input_pdf, output_image, last_name, first_name, middle_initial, topic):
    """ Generate an AI image based on summarized document abstracts. """

    # **Ensure Output Directory Exists**
    output_dir = os.path.dirname(output_image)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # **Step 1: Extract & Summarize Abstracts in Batches**

    abstracts = extract_abstracts_from_pdf(input_pdf, last_name, first_name, middle_initial)
    for i, abstract in enumerate(abstracts, 1):
        print(f"Abstract {i}:\n{abstract}\n")

    # if len(abstracts.split()) < 50:
    #     print("⚠️ Too little text to generate an image. Using topic only.")
    #     abstracts = topic  # Fallback to topic-based generation

    if not abstracts:
        print("⚠️ No abstracts found, using topic only.")
        abstracts = topic  # Fallback

    summarized_text = summarize_in_batches(abstracts, batch_size=3)  # Summarize in batches
    if not summarized_text:
        summarized_text = topic  # Fallback if GPT fails

    # **Step 2: Generate a Structured Prompt**
    structured_prompt = f"""
    Create a collectible-card-game-styled card for {first_name} {last_name}, inspired by their research in {topic}. 
    The card should represent {first_name} {last_name}'s most important contributions.
    Do not include any text on the card.
    Focus on these research themes: {summarized_text}.
    """

    # **Step 3: Generate AI Image with DALL·E**
    try:
        response = client.images.generate(
            model="dall-e-3",
            prompt=structured_prompt,
            size="1024x1024"
        )

        # **Step 4: Save the AI-Generated Image**
        if response and "data" in response:
            image_data = base64.b64decode(response["data"][0]["b64_json"])
            with open(output_image, "wb") as f:
                f.write(image_data)
            print(f"✅ AI-generated image saved as: {output_image}")
            return output_image
        else:
            print("⚠️ No image generated.")
            return None
    except Exception as e:
        print(f"⚠️ Error generating AI image: {e}")
        return None



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

def chat_with_gpt(input_pdf, output_pdf, author, topic):
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
    The search was used to identify articles by the following Author: {author}.
    
    Please summarize the main research themes from this author's work, focusing on the topic of {topic}, since there may be false positives in the search results.

    The extracted text from the PubMed abstracts is as follows:
    ---------------------------
    {pdf_text}
    ---------------------------

    Format the response with:
    - A **Overall Research Themes** by this author section.
    - A **Progress Made in each Field** section.
    - A **section on gaps in knowledge or areas for future research**.
    - If applicable, include **bullet points or numbered lists** for clarity.
    - Ensure clear paragraph separation and professional formatting.
    - Please be clear which articles the findings are based on (i.e., cite the source).

    Begin your response now:
    """

    # **Step 3: Call GPT to Get the Response**
    try:
        stream = client.chat.completions.create(
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


def search_pubmed(last_name, first_name, middle_initial, start_date, end_date, max_results, output_pdf_file1, output_image, topic):
    """ Fetch PubMed articles and save them all into one PDF file, including authors """
    Entrez.email = "your_email@example.com"  # Replace with your email

    first_initial = first_name[0]
    author_query = f'"{last_name} {first_initial}{middle_initial}"[Author]'

    query = f'({author_query}) AND ("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication])'

    # author_query = generate_pubmed_author_query(last_name, first_name, middle_initial)
    # query = f'({author_query}) AND ("{start_date}"[Date - Publication] : "{end_date}"[Date - Publication])'

    print(f"Query: {query}")
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
            # generate_document_image(output_pdf_file1, output_image, last_name, first_name, middle_initial, topic)
        else:
            print("⚠️ No search results found. No PDF generated.")
    else:
        print("No papers found for the specified criteria.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Upload files or directories to S3')
    parser.add_argument('--first', type=str, required=True, help='Author first name to search')
    parser.add_argument('--last', type=str, required=True, help='Last name of author to search')
    parser.add_argument('--middle', type=str, required=True, help='Middle initial of author')
    parser.add_argument('--output', type=str, required=True, help='output first PDF file with results')
    parser.add_argument('--outputImage', type=str, required=True, help='output image')
    parser.add_argument('--topic', type=str, required=True, help='topic of research focus')
    args = parser.parse_args()

    output_pdf_file = args.output
    output_image = args.outputImage
    topic = args.topic
    first = args.first
    last = args.last
    middle = args.middle
    start_date = (datetime.now() - timedelta(days=100000)).strftime("%Y/%m/%d")
    end_date = datetime.now().strftime("%Y/%m/%d")
    max_results = 100  # Maximum number of results to fetch

    search_pubmed(last, first, middle, start_date, end_date, max_results, output_pdf_file, output_image, topic)




