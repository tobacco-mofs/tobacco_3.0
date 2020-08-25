from selenium.webdriver import Chrome
import time
import pandas as pd

#SETTINGS
n_topologies = 3001
webdriver = r"/path/to/chromedriver.exe" #change me!
#^Download from: https://chromedriver.chromium.org/

driver = Chrome(executable_path=webdriver)
topologies = []
df = pd.DataFrame([])
url_tags = ['nets','layers']
print('Scraping topology names..')

for url_tag in url_tags:
    url = "http://rcsr.anu.edu.au/"+url_tag+"#details"
    driver.get(url)
    time.sleep(10)
    driver.find_elements_by_xpath('//*[@id="react-main"]/div/div[2]/ul/li[2]/div/ul/li[1]/a')[0].click()
    for i in range(1,n_topologies+1):
        topologies.append(driver.find_elements_by_xpath('//*[@id="react-main"]/div/div[2]/ul/li[2]/div/div/ul/li['+str(i)+']/a')[0].text)

    print('Scraping topology data...')
    for i, top in enumerate(topologies):

        url = 'http://rcsr.anu.edu.au/nets/'+top
        driver.get(url)
        if i == 0:
            time.sleep(10)
        else:
            time.sleep(1)

        spacegroup = driver.find_elements_by_xpath('//*[@id="react-main"]/div/div[2]/table[1]/tbody/tr/td[2]')[0].text

        cellparams = []
        for j in range(1,7):
            cellparams.append(float(driver.find_elements_by_xpath('//*[@id="react-main"]/div/div[2]/table[2]/tbody/tr/td['+str(j)+']')[0].text))

        n_vertices = int(driver.find_elements_by_xpath('//*[@id="react-main"]/div/div[2]/div[2]/p/span/span[2]')[0].text)
        all_vert_info = []
        for j in range(1,n_vertices+1):
            vert_info = []
            vert_info.append(int(driver.find_elements_by_xpath('//*[@id="react-main"]/div/div[2]/div[2]/table[1]/tbody/tr['+str(j)+']/td['+str(2)+']')[0].text))
            for k in range(3,6):
                vert_info.append(float(driver.find_elements_by_xpath('//*[@id="react-main"]/div/div[2]/div[2]/table[1]/tbody/tr['+str(j)+']/td['+str(k)+']')[0].text))
            all_vert_info.append(vert_info)

        n_edges = int(driver.find_elements_by_xpath('//*[@id="react-main"]/div/div[2]/div[3]/p/span/span[2]')[0].text)
        all_edge_info = []
        for j in range(1,n_edges+1):
            edge_info = []
            for k in range(2,5):
                edge_info.append(float(driver.find_elements_by_xpath('//*[@id="react-main"]/div/div[2]/div[3]/table/tbody/tr['+str(j)+']/td['+str(k)+']')[0].text))
            all_edge_info.append(edge_info)

        df = df.append(pd.DataFrame({'topology':top,'spacegroup':spacegroup,'cellparams':str(cellparams),'vertices':str(all_vert_info),'edges':str(all_edge_info)},index=[0]),ignore_index=True)
        if url_tag == 'nets':
            df.to_csv('rcsr_3D.csv')
        elif url_tag == 'layers':
            df.to_csv('rcsr_2D.csv')

driver.close()
print('Done!')

