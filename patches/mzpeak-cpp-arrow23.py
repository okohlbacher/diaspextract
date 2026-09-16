# Arrow 23 (OpenMS contrib) has only the out-parameter ReadTable; the fork uses the Arrow-24 Result form
# at ONE site. Version-guarded so the source stays correct for Arrow >= 24.
import sys
p = sys.argv[1] if len(sys.argv) > 1 else "mzpeak-cpp/src/util/metadata_model.cpp"
s = open(p).read()
if "ARROW_VERSION_MAJOR < 24" in s: print("already applied"); raise SystemExit
old = """  arrow::Result<std::shared_ptr<arrow::Table>> result =
      metadata.reader().ReadTable();
  if (!result.ok()) {
    throw ParquetError("read metadata table: " + result.status().ToString());
  }
  return std::move(result).ValueOrDie();
"""
new = """#if ARROW_VERSION_MAJOR < 24   // DIAspeXtract cluster build: contrib Arrow 23 has only the out-parameter form
  std::shared_ptr<arrow::Table> table;
  arrow::Status st = metadata.reader().ReadTable(&table);
  if (!st.ok()) {
    throw ParquetError("read metadata table: " + st.ToString());
  }
  return table;
#else
  arrow::Result<std::shared_ptr<arrow::Table>> result =
      metadata.reader().ReadTable();
  if (!result.ok()) {
    throw ParquetError("read metadata table: " + result.status().ToString());
  }
  return std::move(result).ValueOrDie();
#endif
"""
assert s.count(old) == 1, "site not found"
s = s.replace(old, new)
if "#include <arrow/util/config.h>" not in s:
    i = s.find("#include"); s = s[:i] + "#include <arrow/util/config.h>   // ARROW_VERSION_MAJOR\n" + s[i:]
open(p, "w").write(s); print("arrow23 compat applied")
