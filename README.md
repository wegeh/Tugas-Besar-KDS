# Deteksi Mutasi Resistensi Obat pada Gen pol HIV-1

## Deskripsi Singkat

Proyek ini adalah aplikasi web yang bertujuan untuk mendeteksi mutasi resistensi obat pada gen _pol_ HIV-1. Aplikasi menerima input berupa sekuens gen _pol_ pasien dalam format FASTA, melakukan penyelarasan (alignment) terhadap sekuens referensi standar (misalnya HXB2), mengidentifikasi mutasi asam amino, dan kemudian memberikan interpretasi terkait resistensi terhadap berbagai obat antiretroviral (ARV) berdasarkan database aturan mutasi yang telah dikurasi.

Aplikasi ini dibangun menggunakan Python dengan Streamlit untuk antarmuka pengguna, SQLAlchemy untuk interaksi database (Supabase PostgreSQL), serta pustaka bioinformatika seperti Biopython dan Parasail untuk pemrosesan sekuens.

## Fitur Utama

* Unggah sekuens pasien dalam format FASTA.
* Translasi sekuens nukleotida ke protein secara otomatis.
* Pemilihan sekuens referensi dari database.
* Pemilihan metode alignment (Global/Lokal) dan pustaka alignment (Parasail/Biopython).
* Deteksi mutasi asam amino berdasarkan perbandingan dengan sekuens referensi.
* Interpretasi resistensi obat berdasarkan database aturan mutasi yang komprehensif.
* Penyimpanan data sampel, referensi, dan mutasi ke database.
* Ekspor hasil analisis dalam format CSV dan PDF.

## Prasyarat

* Python 3.8 atau lebih tinggi
* Pip (Python package installer)
* Git (opsional, untuk kloning repositori)
* Akses ke database Supabase PostgreSQL (atau PostgreSQL lokal jika Anda mengadaptasi koneksi)

## Pengaturan Proyek

1.  **Kloning Repositori (jika ada):**
    ```bash
    git clone https://github.com/wegeh/Tugas-Besar-KDS.git
    cd Tugas-Besar-KDS
    ```

2.  **Buat dan Aktifkan Virtual Environment:**
    Sangat direkomendasikan untuk menggunakan virtual environment agar dependensi proyek terisolasi.

    * **Windows:**
        ```bash
        python -m venv venv
        .\venv\Scripts\activate
        ```
    * **macOS/Linux:**
        ```bash
        python3 -m venv venv
        source venv/bin/activate
        ```

3.  **Instal Dependensi:**
    ```bash
    pip install -r requirements.txt
    ```

4.  **Konfigurasi Environment Variables:**
    Aplikasi ini memerlukan koneksi ke database Supabase. Salin file `.env.example` menjadi `.env` dan isi detail koneksi database Anda.

    * Salin file .env.example :
        * **Windows:**
            ```bash
            copy .env.example .env
            ```
        * **macOS/Linux:**
            ```bash
            cp .env.example .env
            ```
    * Buka file `.env` yang baru dibuat dan ganti nilai `SUPABASE_DB_URL` dengan string koneksi Supabase PostgreSQL Anda yang valid.

## Menjalankan Aplikasi

1.  **Pastikan Virtual Environment Aktif.**

2.  **Jalankan Aplikasi Streamlit:**
    File utama aplikasi adalah `app.py`.
    ```bash
    streamlit run app.py
    ```
    Aplikasi akan terbuka secara otomatis di browser web default Anda (biasanya di alamat `http://localhost:8501`).

## Struktur Proyek
```
Tugas-Besar-KDS/
├── .venv/                     
├── .env                       
├── .env.example               
├── app.py                     
├── database.py              
├── detect.py             
├── align.py                   
├── translate.py             
├── report.py               
├── requirements.txt           
└── README.md                  
```

## Anggota Kelompok

Berikut adalah anggota kelompok yang berkontribusi dalam proyek ini:

| NIM      | Nama                    |
| :------- | :---------------------- |
| 13522021 | Filbert                 |
| 13522055 | Benardo                 |
| 13522113 | William Glory Henderson |