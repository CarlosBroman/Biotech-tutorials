# Import required packages

import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager

# Install webdriver
driver = webdriver.Chrome(ChromeDriverManager().install())

# Navigate to the webpage and accept cookies.
page_url = "https://pro.boehringer-ingelheim.com/inoncology/our-pipeline"
driver.get(page_url)

WebDriverWait(driver,10).until(EC.element_to_be_clickable((By.XPATH, '//button[text()="Accept All"]'))).click()
WebDriverWait(driver,10).until(EC.element_to_be_clickable((By.XPATH, '//input[@value="I am a Healthcare Professional outside the USA and UK"]'))).click()

page_url = "https://pro.boehringer-ingelheim.com/inoncology/our-pipeline"
driver.get(page_url)

# Extract a list of all the links of the therapies, with names, title and target    
therapy_categories = driver.find_elements(By.CSS_SELECTOR, ".gds-card--visibility-both")
therapies = []

for category in therapy_categories:
    therapy_url = category.find_element(By.CSS_SELECTOR, ".gds-link--default").get_attribute("href")
    if therapy_url.startswith("https://pro.boehringer-ingelheim.com/inoncology/our-pipeline/"):
        therapy_title = category.find_element(By.CSS_SELECTOR, ".gds-card__subtitle").text
        therapy_name = category.find_element(By.CSS_SELECTOR, "h3").text
        therapy_target = category.find_element(By.CLASS_NAME, "cardp").text
        therapies.append({"Therapy": therapy_name, "URL": therapy_url, "Target": therapy_title, "Focus Area": therapy_target})


# Load each therapy page and extract information
for therapy in therapies:
    if "URL" in therapy:
        
        # Get descriptions
        therapy_url = therapy["URL"]
        driver.get(therapy_url)
        therapy_description = driver.find_element(By.XPATH, '//*[@id="main-content"]/div/div/article/div/div[1]/div/div/div/div[2]/div[1]/div/div/p[2]').text

        # Find all the tables on the page
        tables = driver.find_elements(By.CLASS_NAME, "gds-table__wrapper")
        for table in tables:
            # Count the number of rows in the table
            rows = table.find_elements(By.TAG_NAME, "tr")
            num_rows = len(rows)

            # Count the phase I clinical trials
            I_is = 0
            for i in range(0, num_rows):
                cells = rows[i].find_elements(By.TAG_NAME, "td")
                if len(cells) > 1 and cells[1].text.strip() == "I":
                    I_is += 1
            
            # Count the phase II clinical trials
            II_is = 0
            for i in range(0, num_rows):
                cells = rows[i].find_elements(By.TAG_NAME, "td")
                if len(cells) > 1 and cells[1].text.strip() == "II":
                    II_is += 1
                    
            # Count the clinical trials that are atleast partly phase III
            III_is = 0
            for i in range(0, num_rows):
                cells = rows[i].find_elements(By.TAG_NAME, "td")
                if len(cells) > 1 and "III" in cells[1].text.strip():
                    III_is += 1


            # Add the number of clinical trials, the stric phase I and the 
            therapy.update({'Description': therapy_description, 'Clinical Trials': num_rows-1, 'Strictly Phase I': I_is, 'Strictly Phase II': II_is, 'Partially Phase III': III_is})

    else:
        print("Warning: Skipping therapy without URL.")
        
df = pd.DataFrame(therapies)

# Save to csv

df.to_csv("IB_therapies.csv", index=False)