from pathlib import Path
import json,openpyxl
from openpyxl.cell import WriteOnlyCell
from openpyxl.styles import Font,PatternFill,Border,Side,Alignment
from openpyxl.utils import get_column_letter
import argparse
parser=argparse.ArgumentParser(description='Write a private XLSX report from prepared cardiovascular views')
parser.add_argument('--input-json',type=Path,required=True)
parser.add_argument('--output-xlsx',type=Path,required=True)
args=parser.parse_args()
d=json.load(args.input_json.open(encoding='utf8'));w=openpyxl.Workbook(write_only=True)
line=Side(style='thin',color='FFD0D8E0');border=Border(left=line,right=line);headerfill=PatternFill('solid',fgColor='1F4E78');headerfont=Font(name='Arial',size=10,bold=True,color='FFFFFF');font=Font(name='Arial',size=10)
for obj in d['sheets']:
 s=w.create_sheet(obj['name']);s.freeze_panes='F2' if len(obj['headers'])>10 else 'A2';s.sheet_view.showGridLines=False
 for i,h in enumerate(obj['headers'],1):
  width=19
  if h in ('Ген','REF','ALT','Хромосома','Позиция','Ранг','GT','DP','GQ','AD'):width=15
  if any(k in h for k in ('Образцы','HGVS','Consequence','Категории','заболевание','Отбор','контекст','источников','Источник варианта')):width=30
  if h=='Оценки и транскрипты источников':width=35
  s.column_dimensions[get_column_letter(i)].width=width
 s.auto_filter.ref=f'A1:{get_column_letter(len(obj["headers"]))}{len(obj["values"])+1}'
 s.row_dimensions[1].height=32
 head=[]
 for v in obj['headers']:
  c=WriteOnlyCell(s,value=v);c.font=headerfont;c.fill=headerfill;c.border=border;c.alignment=Alignment(vertical='center',wrap_text=True);head.append(c)
 s.append(head)
 for row in obj['values']:
  cells=[]
  for i,v in enumerate(row):
   c=WriteOnlyCell(s,value=v);c.font=font;c.border=border
   if obj['headers'][i] in ('Доля из 4','Доля ALT в AD') and v is not None:c.number_format='0%'
   if obj['headers'][i] in ('gnomAD AF max','gnomAD4 AF','gnomAD MT AF hom','gnomAD MT AF het','Helix MT AF hom','Helix MT AF het') and v is not None:c.number_format='0.00E+00'
   cells.append(c)
  s.append(cells)
 print(obj['name'],len(obj['values']),flush=True)
out=args.output_xlsx;out.parent.mkdir(parents=True,exist_ok=True);w.save(out);print('SAVED',out,flush=True)
