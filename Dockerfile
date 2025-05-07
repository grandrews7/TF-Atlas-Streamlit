FROM python:3.13-slim

COPY requirements.txt .
COPY app.py .

RUN pip install --no-cache-dir -r requirements.txt

EXPOSE 6969

CMD ["streamlit", "run", "app.py", "--server.port", "6969"]